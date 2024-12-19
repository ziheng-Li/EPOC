import PALEOboxes as PB
import PALEOmodel
# import PALEOreactions
import PALEOcopse
import Sundials
import NCDatasets
import Plots

import Infiltrator

"""
    Grid_output_P_O2(
        model, grid_P_norm, grid_pO2PAL; 
        filenameroot, set_P_function, set_O2_function, [tspan, initial_state_use_last] 
    ) -> output_steadystate::PALEOmodel.OutputMemory

Run the supplied model to steady-state over a grid of P, O2 and accumulate the output into `output`

`set_P_function` and `set_O2_function` should be functions that take a single argument (`P_norm` or `pO2PAL`)
and update the model configuration to produce these steady-state values

Returns output_steadystate, a PALEOmodel.OutputMemory with the last timestep (steady-state) for each grid point

Saves 
    - filenameroot*".nc": output_steadystate, last timestep (steady-state) for each grid point
    - filenameroot*"_cube.nc": output_steadystate converted to data-cube by `Grid_cube_P_O2`
"""
function Grid_output_P_O2(
    model::PB.Model,
    grid_P_norm::Vector{Float64},
    grid_pO2PAL::Vector{Float64};
    set_P_function,
    set_O2_function,
    tspan = (0.0, 1e8),
    filenameroot::String,
    initial_state_use_last=false, # true to reuse last grid point as new initial state
)
    
    initial_state, modeldata = PALEOmodel.initialize!(model)

    # create empty output file
    grid_output = PALEOmodel.OutputWriters.OutputMemory()
    nrecords = length(grid_P_norm)*length(grid_pO2PAL)
    PALEOmodel.OutputWriters.initialize!(
        grid_output, model, modeldata, nrecords;
        coords_record=:grid_index,
        coords_record_units="",
    )
   
    #####################################################
    jac_ad=:ForwardDiffSparse # jac_ad=:NoJacobian
    @info "integrate:  ODEProblem using Jacobian $jac_ad"
    f = PALEOmodel.ODE.ODEfunction(
        model, modeldata;
        jac_ad=jac_ad,
        jac_ad_t_sparsity=first(tspan),
        initial_state=initial_state,
    )
    prob = SciMLBase.ODEProblem(f, initial_state, tspan)
    alg = Sundials.CVODE_BDF(linear_solver=:KLU) # alg = Sundials.CVODE_BDF()
    solvekwargs = (reltol=1e-5, saveat=tspan, )
   
    #####################################################
    # Grid for P_norm and pO2PAL steady-state solutions, all the data are saved in output
    #####################################################
    grid_index = 0
    for i in eachindex(grid_pO2PAL)
        pO2PAL = grid_pO2PAL[i]
        set_O2_function(pO2PAL)

        for j in eachindex(grid_P_norm)
            grid_index += 1

            P_norm =grid_P_norm[j]

            println("i=$i (pO2PAL=$pO2PAL) j=$j (P_norm=$P_norm), grid_index=$grid_index")

            set_P_function(P_norm)
            
            # reinitialize model with new pO2PAL, P
            PB.dispatch_setup(model, :setup, modeldata)
            PB.dispatch_setup(model, :initial_value, modeldata)
                       
            @info lpad("", 80, "=")
            @info "integrate: calling solve..."
    
            @time sol = SciMLBase.solve(prob, alg; solvekwargs...);
            PALEOmodel.ODE.print_sol_stats(sol)
            if last(sol.t) < last(tspan)
                @error "run failed at time $(last(sol.t))"
            end
            # @Infiltrator.infiltrate

            PALEOmodel.OutputWriters.add_record!(grid_output, model, modeldata, grid_index)            
            
            if initial_state_use_last
                final_state = PALEOmodel.get_statevar(modeldata.solver_view_all)
                SciMLBase.remake(prob, u0=final_state)
            end
        end
    end

    PALEOmodel.OutputWriters.save_netcdf(grid_output, filenameroot*".nc")

    n_P_norm=length(grid_P_norm)
    n_O2_norm=length(grid_pO2PAL)
    Grid_cube_P_O2(filenameroot, n_P_norm, n_O2_norm)

    return grid_output
end

"""
    Grid_cube_P_O2(filenameroot, n_P_norm, n_O2_norm)

Read `filenameroot*".nc"`, a netcdf with PALEO output over a grid of P, O2
records should be ordered with P changing most rapidly
    
Create an output file `filenameroot*"_cube.nc"` with netcdf dimensions for the P, O2 grid
instead of a single dimension for records,
and map records into this 2D cube
Add coordinates for the P, O2 values (as read from the input file)
"""
function Grid_cube_P_O2(
    filenameroot::String,
    n_P_norm::Int64,    # size of grid - must match input file!
    n_O2_norm::Int64,
)
    input_file = filenameroot*".nc"
    output_file = filenameroot*"_cube.nc"

    n_records = n_P_norm*n_O2_norm

    # data is stored with 
    # grid_index = (i-1)*n_P_norm + j
    # j P_norm varies most rapidly
    # i pO2PAL 

    nc_input =  NCDatasets.NCDataset(input_file, "r")
    nc_output = NCDatasets.NCDataset(output_file, "c")

    try
        # Add coordinates for the P_total, O2 grid
        # NB: there is a value for each record in the input file,
        # we need to "collapse" this to coordinate variables for the 2D grid
        # no error checking that these values actually do correspond to a grid!

        # get O2 grid
        local atm_group = nc_input.group["atm"]
        local O2 = reshape(Array(atm_group["O2"]), (n_P_norm, n_O2_norm))
        local O2_grid = O2[1, :]
        @info "O2_grid: $O2_grid"
        @info "check  : $(O2[end, :])"

        # get P grid
        local ocean_group = nc_input.group["ocean"]
        local P_total = reshape(Array(ocean_group["P_total"]), (n_P_norm, n_O2_norm))
        local P_total_grid = P_total[:, end] # use last, not first record as the high pO2 values are probably more reliable
        @info "P_total_grid: $P_total_grid"
        @info "check  : $(P_total[:, 1])"

        # @Infiltrator.infiltrate

        for (groupname, input_group) in nc_input.group
            @info "groupname $groupname"
            # create group and transfer attributes
            attrib = input_group.attrib
            output_group = NCDatasets.defGroup(nc_output, groupname; attrib=attrib)

            # check we have the correct number of records
            grp_records = input_group.dim["grid_index"]
            grp_records == n_records || error("group records $grp_records != n_records $n_records")

            # transfer dimensions
            for (dimname, dim) in input_group.dim
                @info "    dimension $dimname = $dim"
                NCDatasets.defDim(output_group, dimname, dim)
            end
            # add dimensions and coordinates for P, O grid
            NCDatasets.defDim(output_group, "P_total_grid", n_P_norm)
            NCDatasets.defVar(output_group, "P_total_grid", P_total_grid, ("P_total_grid",))
            NCDatasets.defDim(output_group, "O2_grid", n_O2_norm)
            NCDatasets.defVar(output_group, "O2_grid", O2_grid, ("O2_grid",))

            # transfer variables, changing dimensions from grid_index -> 2D grid
            for (varname, var) in input_group            
                vattrib = Dict(var.attrib)
                # remove netcdf-internal attributes
                filter!(kv -> !(kv.first in ("add_offset", "scale_factor", "_FillValue")), vattrib)
                vdimnames = NCDatasets.dimnames(var)
                vdata = Array(var)
                vsize = size(vdata)
                @info "    varname $varname vdimnames $vdimnames eltype $(eltype(vdata)) vsize $vsize"
                if !isempty(vdimnames) && last(vdimnames) == "grid_index"
                    vdimnames = (vdimnames[1:end-1]..., "P_total_grid", "O2_grid")
                    vsize = (vsize[1:end-1]..., n_P_norm, n_O2_norm)
                    @info "    -> vdimnames $vdimnames vsize $vsize"
                    vdata = reshape(vdata, vsize)
                end
                NCDatasets.defVar(output_group, varname, vdata, vdimnames; attrib)
            end
            
        end
    finally
        close(nc_input)
        close(nc_output)
    end
end

using Plots

"""
    plot_NCDataset_P_O2(output_file; P_total_modern = 3.4e15, P_burial_levels=[2e10, 4e10, 8e10])

Read `output_file` (an offline _cube.nc file produced by `Grid_cube_P_O2`, 
plot P burial flux fluxOceanBurial.flux_total_P as a function of P_norm vs pO2PAL:
- as a heatmap
- as contours at supplied P_burial_levels (these are the critical manifolds with constant P weathering input)
plot Corg burial flux fluxOceanBurial.flux_total_Corg as a function of P_norm vs pO2PAL:
- as a heatmap
- as a line plot vs P_norm (using values at lowest and highest pO2PAL)
"""
function plot_NCDataset_P_O2(
    output_file::String;
    P_total_modern = 3.4e15, # mol P for P_norm=1.0
    P_burial_levels=[2e10, 4e10, 8e10, 12e10], # mol P yr-1
)
    gr(size=(1200, 900))

    NCDatasets.NCDataset(output_file, "r") do ds
        
        O2PAL = 0.21*PB.Constants.k_moles1atm

        group_fluxOceanBurial = ds.group["fluxOceanBurial"]
        P_total_grid = Array(group_fluxOceanBurial["P_total_grid"])
        O2_grid = Array(group_fluxOceanBurial["O2_grid"])

        flux_total_P = Array(group_fluxOceanBurial["flux_total_P"])
        p1 = heatmap(P_total_grid./P_total_modern, O2_grid./O2PAL, flux_total_P'; 
            title="P burial", xlabel="P_total / modern", ylabel="pO2 (PAL)")
        p2 = contour(P_total_grid./P_total_modern, O2_grid./O2PAL, flux_total_P';
            title="P burial $P_burial_levels", xlabel="P_total / modern", ylabel="pO2 (PAL)", levels=P_burial_levels)

        flux_total_Corg = Array(group_fluxOceanBurial["flux_total_Corg"])
        p3 = heatmap(P_total_grid./P_total_modern, O2_grid./O2PAL, flux_total_Corg'; 
            title="Corg burial", xlabel="P_total / modern", ylabel="pO2 (PAL)")
        p4 = plot(P_total_grid./P_total_modern, flux_total_Corg[:, 1],
            title="Corg burial", xlabel="P_total / modern", ylabel="Corg burial (mol yr-1)", label="O2 1")
        
        for i in 1:length(O2_grid)
            plot!(p4, P_total_grid./P_total_modern, flux_total_Corg[:, i],
                title="Corg burial", xlabel="P_total / modern", ylabel="Corg burial (mol yr-1)", label=false) 
        end 

        l = @layout[a b; c d]
        display(plot(p1, p2, p3, p4, layout = l))
    end
end

