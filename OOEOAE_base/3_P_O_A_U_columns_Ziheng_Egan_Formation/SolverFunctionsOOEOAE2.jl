module SolverFunctionsOOEOAE2

import PALEOboxes as PB
import PALEOmodel
import Printf
import Plots
import ColorTypes
import NLsolve
import Interpolations

include("Isoline.jl")
include("ooeoae_expts.jl")

import Infiltrator

"""
    POAnormDeriv(model, modeldata; kwargs...) -> poand::SubsetNormDeriv

See SubsetNormDeriv, calculates normalized derivative and Jacobian for ocean.P, atmocean.O, atmocean.A
holding other variables fixed.

Returns SubsetNormDeriv(model, modeldata,["ocean.P", "atmocean.O", "atmocean.A"]; kwargs...)
"""
function POAnormDeriv(
    model, modeldata; 
    fixed_values=PALEOmodel.get_statevar(modeldata.solver_view_all),
    fixed_t=0.0,
    norm_values=PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
    include_jacobian=false
)
    return SubsetNormDeriv(model, modeldata, ["ocean.P", "atmocean.O", "atmocean.A"]; fixed_values, fixed_t, norm_values, include_jacobian)
end

"""
    POAstatenormDeriv(model, modeldata; kwargs...) -> poand::SubsetNormDeriv

See SubsetNormDeriv, calculates normalized derivative and Jacobian for ocean.P_solve, atmocean.O_solve, atmocean.A_solve
holding other variables fixed.

For those configuration that set reservoirs state_norm = true. ocean.P_solve replaced ocean.P.

Returns SubsetNormDeriv(model, modeldata,["ocean.P", "atmocean.O", "atmocean.A"]; kwargs...)
"""
function POAstatenormDeriv(
    model, modeldata; 
    fixed_values=PALEOmodel.get_statevar(modeldata.solver_view_all),
    fixed_t=0.0,
    norm_values=ones(length(fixed_values)), # ocean.P_solve already being normalized
    include_jacobian=false
)
    return SubsetNormDeriv(model, modeldata, ["ocean.P_solve", "atmocean.O_solve", "atmocean.A_solve"]; fixed_values, fixed_t, norm_values, include_jacobian)
end

"""
    POnormDeriv(model, modeldata; kwargs...) -> poand::SubsetNormDeriv

See SubsetNormDeriv, calculates normalized derivative and Jacobian for ocean.P, atmocean.O,
holding other variables fixed.

Returns SubsetNormDeriv(model, modeldata, ["ocean.P", "atmocean.O"]; kwargs...)
"""
function POnormDeriv(
    model, modeldata; 
    fixed_values=PALEOmodel.get_statevar(modeldata.solver_view_all),
    fixed_t=0.0,
    norm_values=PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
    include_jacobian=false
)
    return SubsetNormDeriv(model, modeldata, ["ocean.P", "atmocean.O"]; fixed_values, fixed_t, norm_values, include_jacobian)
end


"""
    SubsetNormDeriv(model, modeldata, subset_varnames; [fixed_values], [fixed_t], [norm_values], [include_jacobian=false]) -> snd::SubsetNormDeriv

Callable struct that calculates normalized time derivative dv_n/dt(v_n) for subset of scalar variables `subset_varnames`,
eg

    poand = SubsetNormDeriv(model, modeldata, ["ocean.P", "atmocean.O", "atmocean.A"])
    poand(P_n, O_n, A_n) -> (dP_n/dt, dO_n/dt, dA_n/dt)

where O_n is the normalized value of O (ie O = O_n * O_norm where O_norm is read from the provided `norm_values`)

If `include_jacobian=true`, also calculates Jacobian for normalized variables eg P_n, O_n, A_n

    jacobian(poand::SubsetNormDeriv, P_n, O_n, A_n) -> 
        [
            d/dP_n(dP_n/dt)  d/dO_n(dP_n/dt)  d/dA_n(dP_n/dt);
            d/dP_n(dO_n/dt)  d/dO_n(dO_n/dt)  d/dA_n(dO_n/dt);
            d/dP_n(dA_n/dt)  d/dO_n(dA_n/dt)  d/dA_n(dA_n/dt);
        ]
"""
mutable struct SubsetNormDeriv{V1, V2}
    modeldata
    fixed_values::Vector{Float64}
    fixed_t::Float64
    norm_values::Vector{Float64}
    subset_varnames::Vector{String}
    subset_indices::Vector{Int}
    u_work::Vector{Float64}
    du_work::Vector{Float64}
    jac_work::Matrix{Float64}
    ode::V1
    jac_ode::V2
end

function SubsetNormDeriv(
    model, modeldata, subset_varnames; 
    fixed_values=PALEOmodel.get_statevar(modeldata.solver_view_all),
    fixed_t=0.0,
    norm_values=PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
    include_jacobian=false,
)
    iszero(PALEOmodel.num_algebraic_constraints(modeldata.solver_view_all)) || error("algebraic constraints not supported")

    n_full = length(fixed_values)
    n_sub = length(subset_varnames)

    # 20240501 the update of the PB involve the P_solve...
    # if you "state_norm"=>true, then there is no ocean.P but ocean.P_solve instead


    # indices in linear Vector of values used by ODE-like solvers and for Jacobian
    # NB: if an isotope variable (ie multiple indices) only the first index is used
    subset_indices = Int[]
    io = IOBuffer() # only displayed if verbose==true
    println(io, "SubsetNormDeriv: variables and indices in state vector:")
    for vname in subset_varnames
        vindices = PB.get_indices(modeldata.solver_view_all.stateexplicit, vname)
        Printf.@printf(io, "    %20s %s\n", vname, string(vindices))
        length(vindices) <= 2 || error("Variable $vname has more than 2 indices - not a scalar!")
        push!(subset_indices, first(vindices))
    end
    @info String(take!(io))

    ode = PALEOmodel.SolverFunctions.ModelODE(modeldata, modeldata.solver_view_all, modeldata.dispatchlists_all, 0)

    if include_jacobian
        # get ODE-like Jacobian, calculates jac_ode(J, u, p, t) where u is a state Vector
        jac_ode, _ = PALEOmodel.JacobianAD.jac_config_ode(
            :ForwardDiff, model, fixed_values, modeldata, 0.0,
        )
    else
        jac_ode = nothing
    end

    @info String(take!(io))  

    return SubsetNormDeriv(
        modeldata,
        copy(fixed_values),
        fixed_t,
        copy(norm_values),
        copy(subset_varnames),
        subset_indices,
        fill(NaN,  n_full),
        fill(NaN, n_full),
        fill(NaN, n_full, n_full),
        ode,
        jac_ode,
    )
end


function (snd::SubsetNormDeriv)(subset_v_n...)
 
    length(subset_v_n) == length(snd.subset_indices) || 
        error("SubsetNormDeriv number of arguments $(length(subset_v_n)) != length of Variable subset $(length(snd.subset_indices))")

    # create full state vector, from combination of fixed_values and supplied normalized v_n
    snd.u_work .= snd.fixed_values
    for (i, v) in enumerate(subset_v_n)
        i_full = snd.subset_indices[i]
        snd.u_work[i_full] = v*snd.norm_values[i_full]
    end

    # evaluate d/dt 
    snd.ode(snd.du_work, snd.u_work, nothing, snd.fixed_t)

    # return normalized subset derivative
    return snd.du_work[snd.subset_indices] ./ snd.norm_values[snd.subset_indices]
    
end

function jacobian(snd::SubsetNormDeriv, subset_v_n...)
    length(subset_v_n) == length(snd.subset_indices) || 
        error("SubsetNormDeriv number of arguments $(length(v_n)) != length of Variable subset $(length(snd.subset_indices))")

    # create full state vector, from combination of fixed_values and supplied normalized v_n
    snd.u_work .= snd.fixed_values
    for (i, v) in enumerate(subset_v_n)
        i_full = snd.subset_indices[i]
        snd.u_work[i_full] = v*snd.norm_values[i_full]
    end

    # evaluate full Jacobian
    snd.jac_ode(snd.jac_work, snd.u_work, nothing, snd.fixed_t)

    # calculate normalized subset of Jacobian
    Jnorm_subset = zeros(3, 3)
    for i in 1:length(subset_v_n)
        i_full = snd.subset_indices[i]
        for j in 1:length(subset_v_n)
            j_full = snd.subset_indices[j]
            # eg
            # d/dP_n(dP_n/dt), d/dO_n(dP_n/dt), d/dA_n(dP_n/dt)
            # d/dP_n(dO_n/dt), d/dO_n(dO_n/dt), d/dA_n(dO_n/dt)
            # d/dP_n(dA_n/dt), d/dO_n(dA_n/dt), d/dA_n(dA_n/dt)
            Jnorm_subset[i, j] = snd.jac_work[i_full, j_full] * snd.norm_values[j_full] / snd.norm_values[i_full]
        end
    end

    return Jnorm_subset
end

"""
    modify_state(snd::SubsetNormDeriv, subset_v_n...; u_full=copy(snd.fixed_values)) -> u

Modify a state vector, changing indices of subset of variables to supplied normalized values
"""
function modify_state(snd::SubsetNormDeriv, subset_v_n...; u_full=copy(snd.fixed_values))

    for (i, v) in enumerate(subset_v_n)
        i_full = snd.subset_indices[i]
        u_full[i_full] = v*snd.norm_values[i_full]
    end

    return u_full
end

"""
    find_time_series(model, initial_state, modeldata; has_A::Bool, [tspan]) -> (paleorun, A_ts, O_ts, P_ts)

Run paleo model and return time series of A, O, P as Vectors

If has_A=true, assumes model includes A and returns A_ts, otherwise returns A_ts=[]
"""
function find_time_series(model, initial_state, modeldata; has_A::Bool, tspan = (0.0, 1e8), reltol=1e-4, dtmax = 1e4)
    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    # tspan = (0.0, 1e7) # yr
    println("integrate, no jacobian")
    @time PALEOmodel.ODE.integrate(
        paleorun, initial_state, modeldata, tspan, 
        solvekwargs=( # https://diffeq.sciml.ai/stable/basics/common_solver_opts/
            reltol=reltol,
            # reltol=1e-5,
            # saveat=1e6,
            dtmax=dtmax,
        )
    )

    if has_A
        A_ts = PB.get_data(paleorun.output, "atmocean.A_norm")
    else
        A_ts = []
    end
    O_ts = PB.get_data(paleorun.output, "atmocean.O_norm")
    P_ts = PB.get_data(paleorun.output, "ocean.P_norm")

    return (paleorun, A_ts, O_ts, P_ts)
end

"""
    find_nullclines_POA(poand::SubsetNormDeriv; kwargs...) -> (;dPdt_surf, dOdt_surf, dAdt_surf, P_grid, O_grid, A_grid)

Find dP/dt = 0, dO/dt = 0, dA/dt = 0 nullclines in region specified by P_lims, O_lims, A_lims

Each surface is returned as a 2D array on a grid size (n_l, A_npts), with each element a Tuple of P, O, A coordinates,
where n_l is the number of points "across" the surface at constant A, and A_npts is the number of A "levels"

O is treated specially and returns a grid size (n_l, 2) with only the first and last A values.

The arrays of Tuples can be unpacked to separate arrays using `tuples_to_coords(surf)`
"""
function find_nullclines_POA(
    poand::SubsetNormDeriv;
    P_lims = (0.0, 5.0),
    O_lims = (0.0, 2.5),
    A_lims = (2.0, 5.0),
    P_npts = 100, # number of points for P grid
    O_npts = 100, # number of points for O grid
    A_npts = 100,  # number of points for A grid
    n_l = 100, # number of points along each line in dP/dt = 0.0 surface
)
    P_grid = collect(range(P_lims..., length=P_npts))
    O_grid = collect(range(O_lims..., length=O_npts))
    # A_grid = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0] # NB: MarchingCubes needs a regular grid ?!
    A_grid = collect(range(A_lims..., length=A_npts))
   
    # A_grid = [2.0, 2.5]

    dPdt_surf = fill((NaN, NaN, NaN), n_l, length(A_grid))
    dAdt_surf = fill((NaN, NaN, NaN), n_l, length(A_grid))
    # O is a special case, just use the first and last A
    dOdt_surf = fill((NaN, NaN, NaN), n_l, 2)
 

    last_dPdt_line_POA_n_l = nothing # special handling for dP/dt = 0 to avoid glitch in coordinates
    for (k, A_grid_val) in enumerate(A_grid)

        @info "i = $k, A_val = $(A_grid[k])"
        coords = fill((NaN, NaN, NaN), length(P_grid), length(O_grid))

        # create P, O, A coords for P, O plane at constant A
        for i in 1:length(P_grid)
            for j in 1:length(O_grid)
                coords[i, j] = (P_grid[i], O_grid[j], A_grid_val)
            end
        end

        # find isoline dP/dt = 0 in surface defined by coords
        # dP/dt has special handling for glitch in coords using last_dPdt_line_POA_n_l
        # dPdt_line_POA = join_lines_nan_sep(Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[1], coords; maxiters=20)) ## , verbose=true))
        # dPdt_line_POA = only(Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[1], coords; maxiters=20)) ## , verbose=true))
        dPdt_line_POA = Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[1], coords; maxiters=20)[end] ## , verbose=true))
        # interpolate isoline to constant number of points, filling with previous values if coordinates jump (to prevent glitches in the surface plot)
        dPdt_line_POA_n_l = Isoline.interpolate_isoline(dPdt_line_POA, n_l, last_dPdt_line_POA_n_l)     
        
        # # this may cause a error
        last_dPdt_line_POA_n_l = dPdt_line_POA_n_l
        # @info "length of dPdt_line_POA = $(length(dPdt_line_POA)); length of dPdt_line_POA_n_l = $(length(dPdt_line_POA_n_l))"
        # if length(dPdt_line_POA_n_l) != n_l
        #     @error "The length of dPdt_line_POA_n_l should equals to $(n_l)"
        # end

        # fill surface from isolines
        dPdt_surf[:, k] .= dPdt_line_POA_n_l

        # find isoline dA/dt = 0 in surface defined by coords
        dAdt_line_POA = only(Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], coords))
        # to avoid reversals, always put lowest O_n end first
        if first(dAdt_line_POA)[2] > last(dAdt_line_POA)[2]
            reverse!(dAdt_line_POA)
        end
        # interpolate isoline to constant length and fill surface
        dAdt_line_POA_n_l = Isoline.interpolate_isoline(dAdt_line_POA, n_l)
        
        # fill surface from isolines
        dAdt_surf[:, k] .= dAdt_line_POA_n_l

        # find isoline dO/dt = 0 in surface defined by coords
        if A_grid_val == first(A_grid)
            O_k = 1
        elseif A_grid_val == last(A_grid)
            O_k = 2
        else
            O_k = 0
        end
        if !iszero(O_k)
            dOdt_line_POA = Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[2], coords)[]
            # interpolate isoline to constant length and fill surface
            dOdt_line_POA_n_l = Isoline.interpolate_isoline(dOdt_line_POA, n_l)
            # fill surface from isolines
            dOdt_surf[:, O_k] .= dOdt_line_POA_n_l
        end

    end

    return (;
        dPdt_surf,
        dOdt_surf,
        dAdt_surf,
        P_grid,
        O_grid,
        A_grid,
    )

end

"""
    find_nullclines_const_O(poand::SubsetNormDeriv, O_val; kwargs...) -> (;dPdt_line_const_O, dAdt_line_const_O, P_grid, A_grid)

Find dP/dt = 0, dA/dt = 0 nullclines in region specified by P_lims, A_lims at constant O_val
"""
function find_nullclines_const_O(
    poand::SubsetNormDeriv, O_val;
    P_lims = (0.0, 5.0),
    A_lims = (1.0, 5.0),
    P_npts = 100, # number of points for P grid
    A_npts = 100,  # number of points for A grid
)
    P_grid = collect(range(P_lims..., length=P_npts))
    A_grid = collect(range(A_lims..., length=A_npts))

    # create P, O, A coords for P, A plane at constant O
    coords = fill((NaN, NaN, NaN), length(P_grid), length(A_grid))    
    for i in 1:length(P_grid)
        for j in 1:length(A_grid)
            coords[i, j] = (P_grid[i], O_val, A_grid[j])
        end
    end

    # find isoline dP/dt = 0 in surface defined by coords
    # may be more than one line if we've gone outside range, so join them separated by NaN
    dPdt_line_const_O =  join_lines_nan_sep(Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[1], coords))
   
    # find isoline dA/dt = 0 in surface defined by coords
    # may be more than one line if we've gone outside range, so join them separated by NaN
    dAdt_line_const_O = join_lines_nan_sep(Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], coords))
     
    return (;
        dPdt_line_const_O,
        dAdt_line_const_O,
        P_grid,
        A_grid,
    )

end

"""
    find_nullclines_PO(pond::SubsetNormDeriv; [P_lims], [O_lims], [P_npts], [O_npts]) -> (;dPdt_line, dOdt_line, P_grid, O_grid)

Find nullclines dP/dt = 0.0 (the critical manifold) and dO/dt = 0.0 for P, O only model in (P, O) phase plane
"""
function find_nullclines_PO(
    pond::SubsetNormDeriv;
    P_lims = (0.0, 5.0),
    O_lims = (0.0, 2.5),
    P_npts = 100, # number of points for P grid
    O_npts = 100, # number of points for O grid
)
    P_grid = collect(range(P_lims..., length=P_npts))
    O_grid = collect(range(O_lims..., length=O_npts))

    # create P, O coords
    coords = fill((NaN, NaN), length(P_grid), length(O_grid))    
    for i in 1:length(P_grid)
        for j in 1:length(O_grid)
            coords[i, j] = (P_grid[i], O_grid[j])
        end
    end

    # find isoline dP/dt = 0 in surface defined by coords
    # may be more than one line if we've gone outside range, so join them separated by NaN
    dPdt_line =  join_lines_nan_sep(Isoline.find_isolines((P_n, O_n)->pond(P_n, O_n)[1], coords))
   
    # find isoline dO/dt = 0 in surface defined by coords
    dOdt_line = only(Isoline.find_isolines((P_n, O_n)->pond(P_n, O_n)[2], coords))
     
    return (;
        dPdt_line,
        dOdt_line,
        P_grid,
        O_grid,
    )

end


"""
    tuples_to_coords(surf::Array{NTuple{M, Float64}) -> surf_P, surf_O, surf_A

Unpack an Array (or Vector) with a Tuple of M coordinates to M separate coordinate arrays
"""
function tuples_to_coords(surf::Array{NTuple{M, Float64}, N}) where {M, N}    
    surf_c = (map(x->x[m], surf) for m in 1:M)
    return collect(surf_c)
end

"""
    split_line_sign_f(line::AbstractVector, f) -> [(f_is_pos, line_segment), ... ]

Given a line (a Vector of coordinates), split into segments based on value of `f(line[i]) >= 0`
"""
function split_line_sign_f(line, f)
    line_segments = []

    _sign(f, p) = f(p...) >= 0

    current_sign = _sign(f, first(line))
    current_segment = eltype(line)[]
                
    for p in line
        p_sign = _sign(f, p)
        if p_sign != current_sign
            push!(line_segments, (current_sign, current_segment))
            current_segment = eltype(line)[]
            current_sign = p_sign
        end
        push!(current_segment, p)
    end
    push!(line_segments, (current_sign, current_segment))

    return line_segments
end

"""
    plot_segments!(plotfunc, plotfuncfirstarg, linesegs, coord_indices::Vector{Int}, (ls_pos, ls_neg), label; plotfunckwargs...)
"""
function plot_segments!(plotfunc, plotfuncfirstarg, linesegs, coord_indices::Vector{Int}, (ls_pos, ls_neg), label; plotfunckwargs...)
    do_label=label
    for (sign_seg, seg) in linesegs
        seg_coords = SolverFunctionsOOEOAE2.tuples_to_coords(seg)[coord_indices]
        plotfunc(plotfuncfirstarg, seg_coords...; linestyle=sign_seg ? ls_pos : ls_neg, label=do_label, plotfunckwargs...)
        do_label = nothing
    end
end

"""
    join_lines_closest_ends(l1::AbstractVector, l2::AbstractVector)

Concatenate lines, joining at closest ends
"""
function join_lines_closest_ends(l1::AbstractVector, l2::AbstractVector)

    min_dist_between_ends = Inf
    ptsops = [
        # points to test, operations on vector of points before joining
        # l1     l2    op(l1)   op(l2)
        (first, first, reverse, identity),
        (first, last,  reverse, reverse),
        (last,  first, identity, identity),
        (last, last,   identity, reverse),
    ]
    opl1_todo, opl2_todo = nothing, nothing
    for (pl1, pl2, opl1, opl2) in ptsops
        d = Isoline._dist_tuple(pl1(l1), pl2(l2))
        if d < min_dist_between_ends
            opl1_todo, opl2_todo = opl1, opl2
            min_dist_between_ends = d
        end
    end

    return vcat(opl1_todo(l1), opl2_todo(l2))
end

"""
    join_lines_nan_sep(ll::AbstractVector) -> l

Join line segments ll to one line l, separated by NaN
"""
function join_lines_nan_sep(ll::AbstractVector)

    # Infiltrator.@infiltrate
    e = first(first(ll)) # first element of the first line
    sep = Tuple((NaN for i in 1:length(e)))
    # @info "typeof sep = $(typeof(sep))"

    l = copy(first(ll))
    for lseg in ll
        lseg == first(ll) && continue # skip first segment
        push!(l, sep)
        append!(l, lseg)
    end

    return l
end


# reset model, modeldata parameter corg_burial_fac
function reset_corg_burial_fac!(model, modeldata, corgbf_value)

    ooeoae_expts(
        model,  [("corg_burial_fac", corgbf_value)]
    )

    PB.dispatch_setup(model, :setup, modeldata)
end

# solve for corg_burial_fac in range corbf_range that gives dO/dt = 0 at dOdtzeropoint
# pond should be a function such that (dP_n/dt, dO_n/dt) = pond(P_n, O_n)
function find_corg_burial_fac!(model, modeldata, pond, corgbf_range, dOdtzeropoint)

    # find value of dO/dt at (P_n, O_n) = dOdtzeropoint, for given corg_burial_fac
    function dOdtatpoint(corgbf_value, p)
        reset_corg_burial_fac!(model, modeldata, corgbf_value)
        return pond(dOdtzeropoint...)[2]
    end

    # solve for corg_burial_fac that gives dO/dt = 0 at dOdtzeropoint
    corgbf_value, niter, converged = SolverFunctionsOOEOAE2.Isoline.ITP_solve(
        dOdtatpoint, corgbf_range; 
        ϵ=1e-6*sum(corgbf_range)/2,
        verbose=true
    )

    converged || error("find_corg_burial_fac! failed dOdtzeropoint $dOdtzeropoint corgbf_range $cobf_range")

    return corgbf_value
end

# solve for corg_burial_fac in range corbf_range that gives dO/dt = 0 at dOdtzeropoint
# pond should be a function such that (dP_n/dt, dO_n/dt, dA_n/dt) = poand(P_n, O_n, A_n)
function find_corg_burial_fac_line!(model, modeldata, poand, corgbf_range, fold_dPdt_line)

    # find value of dO/dt at (P_n, O_n, A_n) = dOdtzeropoint, for given corg_burial_fac
    function dOdtatdAdt0(corgbf_value, p)
        reset_corg_burial_fac!(model, modeldata, corgbf_value)

        # find where dA/dt=0 nullcline crosses fold line 
        dOdtzeropoint = only(SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], fold_dPdt_line; verbose=true))
        return poand(dOdtzeropoint...)[2]
    end

    # TODO check dAdt really cross the fold -> check is it brackted the root

    # solve for corg_burial_fac that gives dO/dt = 0 at dOdtzeropoint
    corgbf_value, niter, converged = SolverFunctionsOOEOAE2.Isoline.ITP_solve(
        dOdtatdAdt0, corgbf_range; 
        ϵ=1e-6*sum(corgbf_range)/2,
        verbose=true
    )

    converged || error("find_corg_burial_fac! failed corgbf_range $cobf_range")

    return corgbf_value
end

"""
    plot_3D_mapping_to_PO_phase

    P, O phase plane at constant A: find dP/dt = 0 for a representative value of A
    (can get this from surfaces we have already calculated)
"""
function plot_3D_mapping_to_PO_phase(
    dPdt_linesegs_const_A,
    dOdt_dPdt_linesegs_POA, 
    # Obal_P, Obal_O, 
    A_val,
    P_ts, O_ts,
    xlims=[0.0, 5.5],
    ylims=[0.0, 2.25];
    reset_plot_size=true,  # this is a bad idea as it screws up subsequent plots !!
)
    reset_plot_size && Plots.gr(size=(400, 400))
    p1 = Plots.plot()
    SolverFunctionsOOEOAE2.plot_segments!(Plots.plot!, p1, dPdt_linesegs_const_A, [1, 2], (:dash, :solid), "critical manifold, A = $(A_val)"; color=:blue);
    # Plots.plot!(p1, Obal_P, Obal_O; color=:red, label=false); # label="O nullcline"
    SolverFunctionsOOEOAE2.plot_segments!(Plots.plot!, p1, dOdt_dPdt_linesegs_POA, [1, 2], (:dash, :solid), "O_nullcline"; color=:red);

    Plots.plot!(p1, P_ts, O_ts; color=:green, label=false); # label="time series"
    Plots.plot!(p1; xlabel="P_norm", ylabel="O (PAL)", xlims=xlims, ylims=ylims, left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)

    return p1
end

"""
    plot_3D_mapping_to_PA_phase

    P, A plane at constant O: find dP/dt = 0, dA/dt = 0 for one value of O
    (need to calculate nullclines in P, A surface at constant O)
"""
function plot_3D_mapping_to_PA_phase(
    dPdt_linesegs_const_O, 
    dAdt_line_const_O,
    dOdt_dPdt_linesegs_POA,
    O_val,
    P_ts, A_ts;
    O_nullcline::Bool=true,
    A_nullcline::Bool=true,
    xlims=[0.0, 5.5],
    ylims=[0.0, 2.25],
    reset_plot_size=true,  # this is a bad idea as it screws up subsequent plots !!
)
    reset_plot_size && Plots.gr(size=(400, 400))
    p1 = Plots.plot()
    SolverFunctionsOOEOAE2.plot_segments!(Plots.plot!, p1, dPdt_linesegs_const_O, [1, 3], (:dash, :solid), "critical manifold, O = $(O_val)"; color=:blue);
    if A_nullcline
        Plots.plot!(p1, map(x->x[1], dAdt_line_const_O), map(x->x[3], dAdt_line_const_O); color=:black, label="A nullcline");     
    end       
    if O_nullcline
        # Plots.plot!(p1, Obal_P, Obal_O; color=:red, label=false); # label="O nullcline"
        SolverFunctionsOOEOAE2.plot_segments!(Plots.plot!, p1, dOdt_dPdt_linesegs_POA, [1, 2], (:dash, :solid), "O_nullcline"; color=:red);
    end
    Plots.plot!(p1, P_ts, A_ts; color=:green, label="time series");
    Plots.plot!(p1; xlabel="P_norm", ylabel="A_norm", xlims=xlims, ylims=ylims, left_margin = 8Plots.mm, bottom_margin = 8Plots.mm);
    
    return p1
end


"""
    plot_3D_mapping_to_OA_phase

    # O, A plane, "projecting" critical manifold ie using P value(s) from critical manifold at each O, A
"""
function plot_3D_mapping_to_OA_phase(
    folds_dPdt_lines_POA, dOdt_dPdt_linesegs_POA,
    dAdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA_end,
    O_ts, A_ts;
    O_nullcline::Bool=true,
    A_nullcline::Bool=true,
    xlims=[0.0, 5.5],
    ylims=[0.0, 2.25],
    linestyle=:solid,
    reset_plot_size=true,  # this is a bad idea as it screws up subsequent plots !!
)
    reset_plot_size && Plots.gr(size=(400, 400))
    p1 = Plots.plot()
    n_foldlines = length(folds_dPdt_lines_POA);

             # NB: may be 0, 1 or 2 fold lines    
            if n_foldlines == 1
                _, f_O, f_A = SolverFunctionsOOEOAE2.tuples_to_coords(folds_dPdt_lines_POA[1])
                Plots.plot!(p1, f_O, f_A, color=:blue, label="fold")
                Plots.plot!(p1, Plots.Shape(f_O, f_A), color=:blue, alpha=0.1, label=nothing)
            elseif n_foldlines > 1 # just use the first 2
                _, f_O, f_A = SolverFunctionsOOEOAE2.tuples_to_coords(folds_dPdt_lines_POA[1])
                Plots.plot!(p1, f_O, f_A, color=:blue, label="fold")
                _, f_O, f_A = SolverFunctionsOOEOAE2.tuples_to_coords(folds_dPdt_lines_POA[2])
                Plots.plot!(p1, f_O, f_A, color=:blue, label=nothing)
                folds_dPdt_joinclosest_POA = SolverFunctionsOOEOAE2.join_lines_closest_ends(folds_dPdt_lines_POA[1:2]...)
                _, f_O, f_A = SolverFunctionsOOEOAE2.tuples_to_coords(folds_dPdt_joinclosest_POA)
                Plots.plot!(p1, Plots.Shape(f_O, f_A), color=:blue, alpha=0.1, label=nothing)
            end;

            if O_nullcline
                SolverFunctionsOOEOAE2.plot_segments!(Plots.plot!, p1, dOdt_dPdt_linesegs_POA, [2, 3], (linestyle, linestyle), "O_nullcline"; color=:red);
            end

            if A_nullcline
                SolverFunctionsOOEOAE2.plot_segments!(Plots.plot!, p1, dAdt_dPdt_linesegs_POA, [2, 3], (:dash, :solid), "A_nullcline_start"; color=:black);
                if isempty(dAdt_dPdt_linesegs_POA_end)
                else
                    SolverFunctionsOOEOAE2.plot_segments!(Plots.plot!, p1, dAdt_dPdt_linesegs_POA_end, [2, 3], (:dash, :solid), "A_nullcline_end"; color=:black);
                end
            end
            
            Plots.plot!(p1, O_ts, A_ts; color=:green, linestyle=linestyle, label="time series");
            Plots.plot!(p1; xlabel="O_norm", ylabel="A_norm", xlims=xlims, xflip=true, ylims=ylims, left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, legend=:bottomleft)

    return p1
end


"""
linear_inter

Do the linear interpolation, input two values and there corresponding T, input the exact value within the range.
Output the exact T of that point

"""
function linear_inter(
    x_left,
    x_right,
    x_point,
    t_left,
    t_right,
)
    x_span = x_right - x_left
    t_span = t_right - t_left

    t_point = t_left + t_span * (x_point-x_left)/x_span

    return t_point
end

"""
check_factor

In the case the new elements might be the multiple of the existing element, they should be as one
The new input periodic mighe
"""
function check_factor(x, existing_x, threshold)
    is_multiple = false
    for i in [0.5, 1, 2, 3]
        if abs((existing_x * i) - x)/existing_x * 100 < threshold
            is_multiple = true
            break
        end
    end
    return is_multiple
end

function check_error(
    vec,
    threshold, # percantage difference
)
    element_counts = Dict()


    for element in vec
        found = false
        
        for (existing_element, count) in element_counts
            if  check_factor(element, existing_element, threshold)
                element_counts[existing_element] += 1
                found = true
                # @info "same element $element"
                break
            end
        end
        
        if !found
            # @info "different element $element"
            element_counts[element] = 1
        end
    end

    return element_counts
end


"""
find_periodic

Default Input: 
    t_ts, P_ts, O_ts, and
    A specified P_norm value

Output:
    the periodic of the time series, in yr
    the index of the start and end_point
"""
function find_periodic(
    P_ts, O_ts, t_ts, 
    eqb_point; # eqb_point or intersection point, 
    Spec_P = -1, # you can also specify a P_norm value, the function find the starting point for you 
    need_double_check = true,
    min_distance = 0.001,
    omit_t_start=0.0, # time at beginning to omit when estimating P value for zero crossing
    include_t_start_crossings=true, # false to also omit zero crossings where t < omit_t_start
    plot=false,
)

    if need_double_check
        @info """

        ########################### Looking for the periodic #################################
     """
    end

    ####### set default output #######
    (element_counts, start_point_index, end_point_index, sign_change, sign_change_P, sign_change_O) = 
        (Dict(), 1, length(P_ts), [],[],[])
    #######

    if length(eqb_point) == 3
        (eqb_P, eqb_O, _) = eqb_point
    elseif length(eqb_point) <= 2
        (eqb_P, eqb_O) = eqb_point
    end

    eqb_dist = ((eqb_P - P_ts[end])^2 + (eqb_O - O_ts[end])^2)^0.5 

    @info "find_periodic: eqb_dist = $eqb_dist  (eqb_P, eqb_O) = ($eqb_P, $eqb_O), last point = ($(P_ts[end]), $(O_ts[end]))"
    if eqb_dist < min_distance
        @info "The P_norm_end equals to the eqb's P_norm, this might be a 'damping oscillation' case, start from the eqb's P_nrom"
        Spec_P = eqb_P
    elseif Spec_P < 0
        P_ts_filter = Float64[]
        for (t, P) in zip(t_ts, P_ts)
            if t > omit_t_start
                push!(P_ts_filter, P)
            end
        end
        @info "No starting point is been specified, start from the mid of P_norm omitting times < $omit_t_start" 
        Spec_P = (maximum(P_ts_filter) + minimum(P_ts_filter))/2
    else
        @info "Specifing a P_norm = \e[34m$(Spec_P)\e[0m value to find the start point"
        if (Spec_P >= maximum(P_ts)) || (Spec_P <= minimum(P_ts))
            @error "The specified P_norm is out of range"
            return (; element_counts, start_point_index, end_point_index, sign_change, sign_change_P, sign_change_O)
        end
    end

    diff_P = P_ts .- Spec_P
    sign_change_P=[]
    for i in 1:(length(P_ts)-1)
        if include_t_start_crossings || t_ts[i] > omit_t_start
            if (diff_P[i]) * (diff_P[i+1]) < 0 # the point where sign change
                push!(sign_change_P, i)
                # @info "$(diff_P[i]), $(diff_P[i+1])"
            end
        end
    end
    
    if length(sign_change_P) <= 1
        @warn "\e[31mThe specified P_norm = $Spec_P is unique!, this might be a 'globally stable' case"
        return (; element_counts, start_point_index, end_point_index, sign_change, sign_change_P, sign_change_O)
    end

    # # this is the periodic consider both P and O
    # sign_change = intersect(sign_change_P, sign_change_O)

    if plot
        Plots.plot(t_ts, P_ts, label="P_norm")
        Plots.plot!(t_ts, O_ts, label="O_norm")
        display(Plots.scatter!(t_ts[sign_change_P], P_ts[sign_change_P], xlims=(0,5e7)))
    end

    # use the sign change position to find the certain T
    t_vector = []; periodic_vector=[]
    t_vector_P = []; t_vector_O=[]; periodic_vector_P=[]

    for i in 1:2:(length(sign_change_P)-1) 
        t_spec_temp = linear_inter(P_ts[sign_change_P[i]], P_ts[sign_change_P[i]+1], Spec_P,
                            t_ts[sign_change_P[i]], t_ts[sign_change_P[i]+1])
        push!(t_vector_P, t_spec_temp)
    end
    ##########################

    # push the certain periodics into the vector
    # for i in 2:length(t_vector)
    #     push!(periodic_vector, t_vector[i] - t_vector[i-1])
    # end
    for i in 2:length(t_vector_P)
        push!(periodic_vector_P, t_vector_P[i] - t_vector_P[i-1])
    end
    
    if !isempty(periodic_vector_P) # periodic is consistent, stable limit cycle
        (element_counts) = check_error(periodic_vector_P, 5) # 5% error
        (start_point_index, end_point_index) = (sign_change_P[1], sign_change_P[3])
    elseif need_double_check
        @warn "Do not find the periodic, start the double check..."
        pass_double_check = false
        vector_P = vcat(collect(range(start=maximum(P_ts)-0.5, stop=maximum(P_ts), length=5)),
                        collect(range(start=ceil(minimum(P_ts)), stop=floor(maximum(P_ts)), length=5)),
                      collect(range(start=minimum(P_ts), stop=minimum(P_ts)+0.5, length=5)),
                      P_ts[end])
        sort!(vector_P)
        for i in eachindex(vector_P)
            (; element_counts, start_point_index, end_point_index) = SolverFunctionsOOEOAE2.find_periodic(P_ts, O_ts, t_ts, eqb_point; Spec_P=vector_P[i], need_double_check=false)
            if length(element_counts) == 1
                # (start_point_index, end_point_index) = (sign_change[1], sign_change[2])
                pass_double_check = true
            end
        end
        
        if pass_double_check
            @info "\e[32mFinish the double check... we found stable limit cycle\e[0m"
        else
            @info "\e[31mFinish the double check... no stable limit cycle\e[0m"
        end
    else 
        @warn "In the double check process, do not find the periodic..."
    end

    if length(element_counts) == 1
        println("\e[34mLimite cycle is stable, Periodic (yr) = $(first(keys(element_counts))), index = e.g. [$start_point_index, $end_point_index])\e[0m")
    else
        println("\e[31mOscillation is unstable, several Periodics (yr) num = $(length(element_counts)), $(element_counts)\e[0m")
    end

    println("######### Finish find periodic... #########")

    return (; element_counts, start_point_index, end_point_index, sign_change, sign_change_P, sign_change_O)
end

function find_periodic_noinfo(
    P_ts, O_ts, t_ts, 
    eqb_point; # eqb_point or intersection point, 
    Spec_P = -1, # you can also specify a P_norm value, the function find the starting point for you 
    need_double_check = false,
    min_distance = 0.001,
    omit_t_start=0.0, # time at beginning to omit when estimating P value for zero crossing
    include_t_start_crossings=true, # false to also omit zero crossings where t < omit_t_start
    plot=false,
)

    if need_double_check
    #     @info """

    #     ########################### Looking for the periodic #################################
    #  """
    end

    ####### set default output #######
    (element_counts, start_point_index, end_point_index, sign_change, sign_change_P, sign_change_O) = 
        (Dict(), 1, length(P_ts), [],[],[])
    #######

    if length(eqb_point) == 3
        (eqb_P, eqb_O, _) = eqb_point
    elseif length(eqb_point) <= 2
        (eqb_P, eqb_O) = eqb_point
    end

    eqb_dist = ((eqb_P - P_ts[end])^2 + (eqb_O - O_ts[end])^2)^0.5 

    # @info "find_periodic: eqb_dist = $eqb_dist  (eqb_P, eqb_O) = ($eqb_P, $eqb_O), last point = ($(P_ts[end]), $(O_ts[end]))"
    if eqb_dist < min_distance
        # @info "The P_norm_end equals to the eqb's P_norm, this might be a 'damping oscillation' case, start from the eqb's P_nrom"
        Spec_P = eqb_P
    elseif Spec_P < 0
        P_ts_filter = Float64[]
        for (t, P) in zip(t_ts, P_ts)
            if t > omit_t_start
                push!(P_ts_filter, P)
            end
        end
        # @info "No starting point is been specified, start from the mid of P_norm omitting times < $omit_t_start" 
        Spec_P = (maximum(P_ts_filter) + minimum(P_ts_filter))/2
    else
        # @info "Specifing a P_norm = \e[34m$(Spec_P)\e[0m value to find the start point"
        if (Spec_P >= maximum(P_ts)) || (Spec_P <= minimum(P_ts))
            @error "The specified P_norm is out of range"
            return (; element_counts, start_point_index, end_point_index, sign_change, sign_change_P, sign_change_O)
        end
    end

    diff_P = P_ts .- Spec_P
    sign_change_P=[]
    for i in 1:(length(P_ts)-1)
        if include_t_start_crossings || t_ts[i] > omit_t_start
            if (diff_P[i]) * (diff_P[i+1]) < 0 # the point where sign change
                push!(sign_change_P, i)
                # @info "$(diff_P[i]), $(diff_P[i+1])"
            end
        end
    end
    
    if length(sign_change_P) <= 1
        # @warn "\e[31mThe specified P_norm = $Spec_P is unique!, this might be a 'globally stable' case"
        return (; element_counts, start_point_index, end_point_index, sign_change, sign_change_P, sign_change_O)
    end

    # # this is the periodic consider both P and O
    # sign_change = intersect(sign_change_P, sign_change_O)

    if plot
        Plots.plot(t_ts, P_ts, label="P_norm")
        Plots.plot!(t_ts, O_ts, label="O_norm")
        display(Plots.scatter!(t_ts[sign_change_P], P_ts[sign_change_P], xlims=(0,5e7)))
    end

    # use the sign change position to find the certain T
    t_vector = []; periodic_vector=[]
    t_vector_P = []; t_vector_O=[]; periodic_vector_P=[]

    for i in 1:2:(length(sign_change_P)-1) 
        t_spec_temp = linear_inter(P_ts[sign_change_P[i]], P_ts[sign_change_P[i]+1], Spec_P,
                            t_ts[sign_change_P[i]], t_ts[sign_change_P[i]+1])
        push!(t_vector_P, t_spec_temp)
    end
    ##########################

    # push the certain periodics into the vector
    # for i in 2:length(t_vector)
    #     push!(periodic_vector, t_vector[i] - t_vector[i-1])
    # end
    for i in 2:length(t_vector_P)
        push!(periodic_vector_P, t_vector_P[i] - t_vector_P[i-1])
    end
    
    if !isempty(periodic_vector_P) # periodic is consistent, stable limit cycle
        (element_counts) = check_error(periodic_vector_P, 5) # 5% error
        (start_point_index, end_point_index) = (sign_change_P[1], sign_change_P[3])
    elseif need_double_check
        # @warn "Do not find the periodic, start the double check..."
        pass_double_check = false
        vector_P = vcat(collect(range(start=maximum(P_ts)-0.5, stop=maximum(P_ts), length=5)),
                        collect(range(start=ceil(minimum(P_ts)), stop=floor(maximum(P_ts)), length=5)),
                      collect(range(start=minimum(P_ts), stop=minimum(P_ts)+0.5, length=5)),
                      P_ts[end])
        sort!(vector_P)
        for i in eachindex(vector_P)
            (; element_counts, start_point_index, end_point_index) = SolverFunctionsOOEOAE2.find_periodic_noinfo(P_ts, O_ts, t_ts, eqb_point; Spec_P=vector_P[i], need_double_check=false)
            if length(element_counts) == 1
                # (start_point_index, end_point_index) = (sign_change[1], sign_change[2])
                pass_double_check = true
            end
        end
        
        # if pass_double_check
        #     @info "\e[32mFinish the double check... we found stable limit cycle\e[0m"
        # else
        #     @info "\e[31mFinish the double check... no stable limit cycle\e[0m"
        # end
    else 
        # @warn "In the double check process, do not find the periodic..."
    end

    if length(element_counts) == 1
        # println("\e[34mLimite cycle is stable, Periodic (yr) = $(first(keys(element_counts))), index = e.g. [$start_point_index, $end_point_index])\e[0m")
    else
        # println("\e[31mOscillation is unstable, several Periodics (yr) num = $(length(element_counts)), $(element_counts)\e[0m")
    end

    # println("######### Finish find periodic... #########")

    return (; element_counts, start_point_index, end_point_index, sign_change, sign_change_P, sign_change_O)
end

"""
    find_periodic_wapper

    The PALEO wapper for find_periodic,
    The initial input is the model output
"""
function find_periodic_wapper(output::PALEOmodel.OutputWriters.OutputMemory)
    t_ts = PB.get_data(output, "global.tforce")
    P_ts = PB.get_data(output, "ocean.P_norm")
    O_ts = PB.get_data(output, "atmocean.O_norm")
    eqb_point = (P_ts[end], O_ts[end])

    return find_periodic(P_ts, O_ts, t_ts, eqb_point)
end

"""
    find_sharpness_corg_burial_fac

    Input: vector_pars, 
    
    Output: sharpness, corg_burial_fac
"""
function find_sharpness_corg_burial_fac(vector_pars)
    sharpness = 1000
    corg_burial_fac = 1000
    CP_ratio = -1

    # we need to update the CP_ratio first
    for expt in vector_pars
        if length(expt) == 3 && expt[1] == "CPsea"
            _, k_oxic, k_anoxic = expt
            CP_ratio = k_anoxic/k_oxic
        end
    end

    for expt in vector_pars
        if length(expt) == 2 && first(expt) == "k_O2_U"
            _, (k_O2_U_min, k_O2_U_max) = expt
            if CP_ratio < 0
                @error "find_sharpness_corg_burial_fac need the CPsea in the vector_pars!"
                return (; sharpness, corg_burial_fac)
            end
            sharpness = find_sharpness((k_O2_U_min, k_O2_U_max), CP_ratio)

        elseif length(expt) ==2 && first(expt) == "corg_burial_fac"
            _, corg_burial_fac = expt
        end
    end

    return (; sharpness, corg_burial_fac)
end

function find_sharpness((k_O2_U_min, k_O2_U_max), CP_ratio)
    # sharpness = k_O2_U_min * CP_ratio - k_O2_U_max
    # SJD normalize to range 0 - 1
    k_O2_U_mean = 0.5*(k_O2_U_min + k_O2_U_max)
    sharpness = (k_O2_U_min * CP_ratio - k_O2_U_max)/(k_O2_U_mean*(CP_ratio-1))

    return sharpness
end

function find_sharpness_PALEOwapper(output::PALEOmodel.OutputWriters.OutputMemory, CP_ratio)
    O2_U = PB.get_data(output, "ocean.O2_U_local")

    sharpness=[]
    for i in 1:length(O2_U)
        push!(sharpness, find_sharpness((O2_U[i][1], O2_U[i][end]), CP_ratio))
    end

    return sharpness
end



"""
    find_intersection_point(model, modeldata, poand, corg_burial_fac, dPdt_surf) -> (P_n, O_n, A_n)

Find intersection point for specified corg_burial_fac
    
Assumes dPdt_surf doesn't change, then recalculates O, A nullclines in that surface and finds intersection

NB: resets model to supplied corg_burial_fac parameter !!!
"""
function find_intersection_point(model, modeldata, poand, corg_burial_fac, dPdt_surf)
    reset_corg_burial_fac!(model, modeldata, corg_burial_fac)
    # Recalculate nullclines in dP/dt=0 surface
    # dO/dt = 0 intersection with dP/dt = 0
    dOdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[2], dPdt_surf))
    # find where dA/dt=0 along line with dO/dt=0 and dP/dt=0
    intersection = only(SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dOdt_dPdt_line_POA; verbose=true))

    return intersection
end


"""
    find_intersection_point_newton(poand) -> (P_n, O_n, A_n)

Find intersection point (stable or unstable eqb point, where d/dt = 0) using Newton's method

Supplied `poand` should provide:

    poand(P_n, O_n, A_n) -> (dP_n/dt, dO_n/dt, dA_n/dt)
    SolverFunctionsOOEOAE2.jacobian(poand, P_n, O_n, A_n) -> jacobian

NB: requires forked version of NLsolve, add this to the PALEOexamples environment using:

    (PALEOexamples) pkg> add https://github.com/PALEOtoolkit/NLsolve.jl#project_region
"""
function find_intersection_point_newton(poand; initialpoa=[1.0, 1.0, 1.0])
    
    # define adaptor functios for NLsolve (requires vector [P_n, O_n, A_n])
    # define a function to calculate f (time derivative)
    f(POA_n) = poand(POA_n...)
    # define a function to calculate Jacobian
    jac(POA_n) = SolverFunctionsOOEOAE2.jacobian(poand, POA_n...)
    
    sol = NLsolve.nlsolve(
        f, jac, initialpoa;
        ftol=1e-12,
        method=:newton,
        linesearch=NLsolve.LineSearches.Static(),
        store_trace=true,
        show_trace=true,
        extended_trace=true,
        iterations=100,
        apply_step! = PALEOmodel.SolverFunctions.StepClampMultAll!(1e-3, 5.0, 1.2)
    )

    # @Infiltrator.infiltrate
    # sol.f_converged == false, sol.iterations == 100


    return (sol.zero...,)
end

"""
    resample_map_grid(x_old_grid, y_old_grid, z_old, x_new_left, x_new_right; kwargs...)
        -> x_new_grid, y_new_grid, z_new

Resample and map a 2D array z_old with grid coordinates x_old_grid, y_old_grid to z_new with
grid coordinates x_new_grid, y_new_grid

The x coordinates are transformed where x_old_grid is assumed to cover the interval (0, 1) which is
linearly transformed into the interval (x_new_left(y), x_new_right(y)), ie

    (0, y) -> (x_new_left(y), y) 
    (1, y) -> (x_new_right(y), y)

and resampled so z_new has n_x_new regularly spaced points (default 200)

The y coordinates are unchanged, and are resampled so z_new has n_y_new regularly spaced points (default 200)
"""
function resample_map_grid(
    x_old_grid, y_old_grid, z_old, x_new_left, x_new_right;
    n_x_new=200, x_new_grid_left=minimum(x_new_left),  x_new_grid_right=maximum(x_new_right),
    n_y_new=200,
)
    size(z_old, 1) == length(y_old_grid) || error("z_old y_old_grid sizes incorrect")
    size(z_old, 2) == length(x_old_grid) || error("z_old x_old_grid sizes incorrect")
    x_new_grid_left < x_new_grid_right || error("x_new_right must be larger than x_new_left")

    x_new_grid = [x_new_grid_left + (j-1)/(n_x_new-1)*(x_new_grid_right-x_new_grid_left) for j in 1:n_x_new]
    y_new_grid = [first(y_old_grid) + (i-1)/(n_y_new-1)*(last(y_old_grid)-first(y_old_grid)) for i in 1:n_y_new]
    y_new_grid[end] = y_old_grid[end] # floating point gotcha:  y_new_grid may be epsilon outside y_old_grid which confuses naive interpolation

    # Infiltrator.@infiltrate

    z_new = zeros(n_y_new, n_x_new)

    # 2D linear interpolation
    # interpolation in x is handled by Interpolations
    # interpolation in y we do ourselves
    z_old_interp_x = [Interpolations.linear_interpolation(x_old_grid, z_old[i, :]) for i in 1:length(y_old_grid)]

    for (i, y_new) in enumerate(y_new_grid)
        # y interpolation
        # find the two rows in z_old that we will linearly interpolate between
        i_old_upper, i_old_lower = -1, -1  # indices of rows
        ai_old_upper, ai_old_lower = NaN, NaN # linear interpolation coeffs
        for io_upper in 2:length(y_old_grid)
            if (y_new - y_old_grid[io_upper-1])*(y_old_grid[io_upper] - y_new) >= 0 # y_new in interval or at endpoint
                i_old_lower, i_old_upper = io_upper-1, io_upper
                y_old_lower, y_old_upper = y_old_grid[i_old_lower], y_old_grid[i_old_upper]
                ai_old_upper = (y_new - y_old_lower)/(y_old_upper - y_old_lower)
                ai_old_lower = 1.0 - ai_old_upper
                break
            end
        end

        # find the new x coordinates corresponding the x_old = (0, 1), for this y_new
        x_new_0 = ai_old_lower*x_new_left[i_old_lower] + ai_old_upper*x_new_left[i_old_upper]
        x_new_1 = ai_old_lower*x_new_right[i_old_lower] + ai_old_upper*x_new_right[i_old_upper]

        # fill x row
        for (j, x_new) in enumerate(x_new_grid)
            x_old = (x_new - x_new_0)/(x_new_1 - x_new_0)

            if (x_old < first(x_old_grid)) || (x_old > last(x_old_grid))
                z_new[i, j] = NaN
            else
                # 2D linear interpolation
                # interpolation in x is handled by Interpolations
                # interpolation in y we do ourselves
                z_new[i, j] = (
                    ai_old_lower*z_old_interp_x[i_old_lower](x_old) +
                    ai_old_upper*z_old_interp_x[i_old_upper](x_old)
                )
            end
        end
    end

    return x_new_grid, y_new_grid, z_new
end

"""
    forcing_white_noise

    generate white noise, norm_distribution return a time series

    The input:  
        step: time interval between two pulse
        dt: duration of each pulse

    The output of this function is prepared for the Reaction CO2pulse, perturb_totals and perturb_times
"""
function forcing_white_noise(
    start::Float64,    # e.g. -600e6
    step::Float64,     # e.g. 1e6, step of the time series, 1e6 per point
    stop::Float64,     # e.g. -300e6
    CO2pulse::Float64; # e.g. 1e14 mol/yr
    dt::Float64=1e2,   # e.g. 1e4
)
    if dt > step
        @error "The duration of the pulse should never greater than between two pulses"
    end

    age = collect(range(start=start,step=step,stop=stop))
    pulse = randn(length(age)) * CO2pulse     

    # age_temp = collect(range(start=start,step=step,stop=stop))
    # pulse_temp = abs.(randn(length(age_temp))) * CO2pulse # positive part of norm_distribution * CO2pulse
    # age=Float64[]
    # pulse=Float64[]

    # for i in eachindex(age_temp)
    #     age = vcat(age, [age_temp[i]-1.0, age_temp[i], age_temp[i]+dt, age_temp[i]+dt+1.0])
    #     pulse = vcat(pulse, [0.0, pulse_temp[i], pulse_temp[i], 0.0])
    # end

    return age, pulse
end 

end # module