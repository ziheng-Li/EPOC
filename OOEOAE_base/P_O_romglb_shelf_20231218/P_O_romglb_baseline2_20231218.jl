import PALEOboxes as PB
import PALEOmodel
# import PALEOreactions
import PALEOcopse
import PALEOocean
import Sundials
import SciMLBase
import DataFrames
using Plots

#############################
# Test shelf Corg burial vs O2 utilization
# 2nd version: use true oxygen utilisation from a run with zero nutrients
#############################

gr(size=(400, 400))
# gr(size=(600, 600))


include("../ooeoae_expts.jl")
include("../ooeoae_plots.jl")
# include("../GLODAPv2.2020/base_GLODAP.jl")

include("../../PALEOreactions/OceanTransportRomanielloShelf.jl")
include("../../PALEOreactions/Burial.jl")
include("../../PALEOreactions/AtmReservoirs.jl")

output_folder_name = "P_O_romglb_baseline_20231218"
# dropbox_output_dir = "/home/sd336/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/P_O_romglb_shelf_20231116" # will need to create this manually
# dropbox_output_dir = "/Users/liziheng/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/P_O_romglb_shelf_20231116"
dropbox_output_dir = joinpath(@__DIR__, "../../figures/P_O_romglb_shelf_20231116") # will need to create this manually

output_figures_dir = joinpath(dropbox_output_dir, output_folder_name)
isdir(output_figures_dir) || mkdir(output_figures_dir)

#####################################################
# Create model
#####################################################

# define a target cumulative Corg burial vs O2_U distribution
# as a linear function of O2_U 
# target_Corg_total = 8e12 # mol Corg yr-1
target_Corg_total = 4.5e12 # mol Corg yr-1

tspan = (0.0, 1e6) # yr

# NB: nutrient input into all surface boxes increases shelf nutrients and will fail if a shelf area is ~0,
# so use oceanhlatgyresurface for nutrient input

expts_table = [
    (
        "burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p23_0p25", ["../transportromglbshelf7.yaml","../P_O_romglb_hlatgyre.yaml"], "model_BioProdPrest_Martin_Ozaki2011", 
        (0.23, 0.25), # O2_U min, max
        target_Corg_total - 0.44e12, # target for shelf adjustment
        [
            ("biopumpCorg_Martin", true, 0.858, 100.0),
            # ("toceanexch", 2 .*vcat([3.0, 3.0, 6.0], fill(3.0, 19))),  # NB: 1st shelf adjusted to be more oxic
            # ("toceanexch2", 2 .* vcat([6.0, 3.0, 3.0], fill(3.0, 19))), 
            ("tshelfexch", fill(0.5, 21)),
            # shelf_Afloor_adj calculated from first run with no area adjustment
            ("shelf_areas", [1.0e6, 1.0e6, 1.0e6, 1.9689093257290117e13, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]),
            # 1st iteration
            # ("shelf_areas", [1.0e6, 1.0e6, 3.937307327229012e13, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]),
            ("VCI_Bergman", 0.124, 5e-3, 6e-3, 112.5*3, 112.5*3, 112.5*3, 4), # oxic P burial 4e10, Corg burial 4.5e12
        ],
    ),
    (
        "burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p15_0p6", ["../transportromglbshelf7.yaml","../P_O_romglb_hlatgyre.yaml"], "model_BioProdPrest_Martin_Ozaki2011", 
        (0.15, 0.6), # O2_U min, max
        target_Corg_total,
        [
            ("biopumpCorg_Martin", true, 0.858, 100.0),
            # ("toceanexch", 2 .*vcat([3.0, 3.0, 6.0], fill(3.0, 19))),  # NB: 1st shelf adjusted to be more oxic
            # ("toceanexch2", 2 .* vcat([6.0, 3.0, 3.0], fill(3.0, 19))),
            ("tshelfexch", fill(0.5, 21)),
            # shelf_Afloor_adj calculated from first run with no area adjustment
            # no area adjustment
            ("shelf_areas", [1.0e6, 1.3744579651577837e12, 3.681906566607758e12, 5.154563044616292e12, 2.136323704396951e12, 1.3667499775004841e12, 1.3181670106043784e12, 1.1809753978705688e12, 9.88421539255121e11, 5.385913978724246e11, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]),
            # 1st iteration
            # ("shelf_areas", [1.0e6, 1.0e6, 1.1643002254645312e13, 5.059211024835465e12, 2.6191164172280103e12, 2.393159126131237e12, 4.05880110600462e12, 1.8948922267373552e12, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]),
            ("VCI_Bergman", 0.124, 5e-3, 6e-3, 112.5*3, 112.5*3, 112.5*3, 4), # oxic P burial 4e10, Corg burial 4.5e12
        ],
    ),
    (
        "burial_shelves6_BioProdPrest_Martin_Ozaki2011_0p05_0p75", ["../transportromglbshelf7.yaml","../P_O_romglb_hlatgyre.yaml"], "model_BioProdPrest_Martin_Ozaki2011", 
        (0.05, 0.75), # O2_U min, max
        target_Corg_total,
        [
            ("biopumpCorg_Martin", true, 0.858, 100.0),
            # ("toceanexch", 2 .*vcat([3.0, 3.0, 6.0], fill(3.0, 19))),  # NB: 1st shelf adjusted to be more oxic
            # ("toceanexch2", 2 .* vcat([6.0, 3.0, 3.0], fill(3.0, 19))),
            ("tshelfexch", fill(0.5, 21)),
            # shelf_Afloor_adj calculated from first run with no area adjustment
            ("shelf_areas", [5.32783275113773e12, 3.244095365400866e12, 2.3669399356764155e12, 3.3136476715390405e12, 1.373350952826603e12, 7.108272527441417e11, 6.329331221523818e11, 6.856142292673103e11, 6.354138466640063e11, 4.7935543009323285e11, 4.196450264168834e11, 3.505760694687173e11, 3.150484359984442e11, 2.5680100070643924e11, 2.1701895211484616e11, 1.6950067992189218e11, 1.2688915840575194e11, 1.0e6, 1.0e6, 1.0e6, 1.0e6]),
            # 1st iteration
            # ("shelf_areas", [1.0e6, 7.859620559812161e12, 7.611402752239954e12, 3.399149843991235e12, 1.515453829446453e12, 1.0522045886743362e12, 2.5407903104201157e12, 1.4004074964471104e12, 8.329845096885391e11, 7.22017404252616e11, 5.157699622075906e11, 4.3122575808480347e11, 5.744156828060865e11, 3.3384208203759924e11, 8.59088682780766e10, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]),
            ("VCI_Bergman", 0.124, 5e-3, 6e-3, 112.5*3, 112.5*3, 112.5*3, 4), # oxic P burial 4e10, Corg burial 4.5e12
        ],
    ),
]

# box indices for ../transportromglbshelf4.yaml
ihlat =     collect(1:13) # (surface box is 1)
igyre =    collect(14:46)
iupw  =     collect(47:79)
# nshelf = 18
nshelf = 21
ilast_shelf = 80 + 2*nshelf-1
ishelfX =   collect(80:ilast_shelf)
ishelfX_surf =   collect(80:2:ilast_shelf)
ishelfX_floor = ishelfX_surf .+ 1

P_O_romglb_Table5 = Dict() # all results, indexed by fileroot
P_O_romglb_shelf_Afloor_adj = Dict() # adjusted shelf areas, indexed by fileroot

p_corg_cumul_summary = plot(;
    xlabel="Accumulated Corg burial (mol yr-1)", ylabel="frac O2 utilization",
    bottom_margin = 1Plots.mm, right_margin = 1Plots.mm,
)

p_shelf_area_summary = plot(;
    xlabel="shelf number", ylabel="Shelf area (m^2)",
    yscale=:log10,
    bottom_margin = 1Plots.mm, right_margin = 1Plots.mm,
)

# the simplest way to use consistent colors is to refer to them by index in the current "cycle"
# https://discourse.julialang.org/t/way-to-get-default-color-order-in-plots-jl/12643/11
coloridx = 1

#####################################################
# df_hlat = Data_filter(["Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 6000.0, "highlat", "oxygen")
# df_gyre = Data_filter(["Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 6000.0, "lowlat", "oxygen")
# df_upw = Data_filter(["Arctic_Ocean" "Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 6000.0, "all", "oxygen", 10)
# df_shelf = Data_filter(["Arctic_Ocean" "Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 200.0, "all", "oxygen")

for (fileroot, yamls, model_name, O2_U_target, target_Corg_adjust, vector_pars) in expts_table
# (fileroot, yamls, model_name, O2_U_target, vector_pars) = only(expts_table)

    local target_O2_U_min, target_O2_U_max = O2_U_target


    #####################################################
    # Create model
    #####################################################
    local model = PB.create_model_from_config(
        joinpath.(@__DIR__, yamls), 
        model_name, 
        modelpars=Dict(),
    )

    ooeoae_expts(
        model, vector_pars
    )
    local initial_state, modeldata = PALEOmodel.initialize!(model)

    local paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    P_O_romglb_Table5[fileroot] = paleorun

    println("integrate, no jacobian")
    @time PALEOmodel.ODE.integrate(
        paleorun, initial_state, modeldata, tspan, 
        solvekwargs=( # https://diffeq.sciml.ai/stable/basics/common_solver_opts/
            reltol=1e-4,
        )
    )
    
    # gr(size=(650, 400))
    local p1 = plot_vs_depth(paleorun.output, "ocean.Abox",                   "Ocean area (m2)", "Depth (m)", NaN, "", true)
    local p2 = plot_vs_depth(paleorun.output, "ocean.bioprod/Prod_Corg",      "Bioprod (mol/yr)", "Depth (m)", NaN, "", true)
    local p3 = plot_vs_depth(paleorun.output, "ocean.remin_Corg",             "Remin Corg (mol/yr)", "Depth (m)", NaN, "", true)
    # local p3 = plot_vs_depth(paleorun.output, "ocean.export_Corg", "Export Corg (mol/yr)", "Depth (m)", NaN, "", true) # export = bioprod
    local p4 = plot_vs_depth(paleorun.output, "fluxOceanBurial.flux_Corg",    "Corg burial (mol/yr)", "Depth (m)", NaN, "", true)
    local l = @layout[a b; c d]
    local p = plot(p1, p2, p3, p4, layout = l)
    display(p)
    savefig(joinpath(output_figures_dir, "$(fileroot)_vs_depth.svg"))

    # plot_compare_GLODAP(paleorun.output, df_hlat, df_gyre, df_upw, df_shelf, "plot_table5/$(fileroot).png")

    #####################################################
    # Create model with no nutrients
    #####################################################
    local model_noP = PB.create_model_from_config(
        joinpath.(@__DIR__, yamls), 
        model_name, 
        modelpars=Dict(
            "P_restore"=>0.0,
        ),
    )

    ooeoae_expts(
        model_noP, vector_pars
    )
    local initial_state, modeldata = PALEOmodel.initialize!(model_noP)

    local paleorun_noP = PALEOmodel.Run(model=model_noP, output = PALEOmodel.OutputWriters.OutputMemory())

    # P_O_romglb_Table5[fileroot] = paleorun

    println("integrate, no jacobian")
    @time PALEOmodel.ODE.integrate(
        paleorun_noP, initial_state, modeldata, tspan, 
        solvekwargs=( # https://diffeq.sciml.ai/stable/basics/common_solver_opts/
            reltol=1e-4,
        )
    )

    #################################
    # Sum plot number of shelves vs. O2_U (single panels)
    #################################

    local df = DataFrame(
        i=1:last(ishelfX_floor), # box index
        ocean_zlower = PB.get_data(paleorun.output, "ocean.zlower")[end],
        oceanfloor_Afloor = PB.get_data(paleorun.output, "oceanfloor.Afloor")[end],
        ocean_O2_conc = PB.get_data(paleorun.output, "ocean.O2_conc")[end],
        ocean_O2_surf_conc = PB.get_data(paleorun.output, "ocean.O2_conc")[end], # placeholder
        ocean_P_conc = PB.get_data(paleorun.output, "ocean.P_conc")[end],
        Corg_b = PB.get_data(paleorun.output, "fluxOceanBurial.flux_Corg")[end],
    )

    local df_noP = DataFrame(
        i=1:last(ishelfX_floor), # box index
        ocean_zlower = PB.get_data(paleorun_noP.output, "ocean.zlower")[end],
        oceanfloor_Afloor = PB.get_data(paleorun_noP.output, "oceanfloor.Afloor")[end],
        ocean_O2_conc = PB.get_data(paleorun_noP.output, "ocean.O2_conc")[end],
        ocean_O2_surf_conc = PB.get_data(paleorun_noP.output, "ocean.O2_conc")[end], # placeholder
        ocean_P_conc = PB.get_data(paleorun_noP.output, "ocean.P_conc")[end],
        Corg_b = PB.get_data(paleorun_noP.output, "fluxOceanBurial.flux_Corg")[end],
    )

    df.ocean_depth = - df.ocean_zlower
    # add [O2] for corresponding surface box and calculate fractional O2 utilisation
    df.ocean_O2_surf_conc[ihlat] .= df.ocean_O2_surf_conc[first(ihlat)] 
    df.ocean_O2_surf_conc[igyre] .= df.ocean_O2_surf_conc[first(igyre)] 
    df.ocean_O2_surf_conc[iupw] .= df.ocean_O2_surf_conc[first(iupw)] 
    df.ocean_O2_surf_conc[ishelfX_floor] .= df.ocean_O2_surf_conc[ishelfX_surf] 
    df.ocean_O2_U = (df_noP.ocean_O2_conc .- df.ocean_O2_conc)./df_noP.ocean_O2_conc

    # sort!(df, :ocean_zlower)

    local df_all_O2_conc = copy(df)
    sort!(df_all_O2_conc, :ocean_O2_conc)
    df_all_O2_conc.Corg_b_cumsum = cumsum(df_all_O2_conc.Corg_b)
    p = plot(
        df_all_O2_conc.Corg_b_cumsum, df_all_O2_conc.ocean_O2_conc;
        xlabel="Accumulated Corg burial (mol yr-1)", ylabel="[O2] (mol m-3)", label="all",
        bottom_margin = 1Plots.mm, right_margin = 1Plots.mm,
    )
    display(p)

    local df_all_O2_U = copy(df)
    sort!(df_all_O2_U, :ocean_O2_U)
    df_all_O2_U.Corg_b_cumsum = cumsum(df_all_O2_U.Corg_b)
    p = plot(
        df_all_O2_U.Corg_b_cumsum, df_all_O2_U.ocean_O2_U;
        xlabel="Accumulated Corg burial (mol yr-1)", ylabel="frac O2 utilization", label="all",
        bottom_margin = 1Plots.mm, right_margin = 1Plots.mm,
    )
    display(p)

    local df_shelves_O2_conc = df[ishelfX_floor, :]
    sort!(df_shelves_O2_conc, :ocean_O2_conc)
    df_shelves_O2_conc.Corg_b_cumsum = cumsum(df_shelves_O2_conc.Corg_b)
    p = plot(
        df_shelves_O2_conc.Corg_b_cumsum, df_shelves_O2_conc.ocean_O2_conc;
        xlabel="Accumulated Corg burial (mol yr-1)", ylabel="[O2] (mol m-3)", label="shelves",
        xlims=(0.0, Inf), ylims=(0.0, Inf),
        bottom_margin = 1Plots.mm, right_margin = 1Plots.mm,
    )
    display(p)

    local df_shelves_O2_U = df[ishelfX_floor, :]
    sort!(df_shelves_O2_U, :ocean_O2_U)
    df_shelves_O2_U.Corg_b_cumsum = cumsum(df_shelves_O2_U.Corg_b)
    p = plot(
        df_shelves_O2_U.Corg_b_cumsum, df_shelves_O2_U.ocean_O2_U;
        xlabel="Accumulated Corg burial (mol yr-1)", ylabel="frac O2 utilization", label="shelves",
        xlims=(0.0, Inf), ylims=(0.0, Inf),
        bottom_margin = 1Plots.mm, right_margin = 1Plots.mm,
    )
    display(p)

    # unadjusted total Corg burial
    Corg_b_all = last(df_all_O2_U.Corg_b_cumsum)
    Corg_b_shelves = last(df_shelves_O2_U.Corg_b_cumsum)
    Corg_b_notshelves = Corg_b_all - Corg_b_shelves

    function target_Corg_cumul(O2_U)  
        O2_U_clamp = clamp(O2_U, target_O2_U_min, target_O2_U_max)
        Corg_cumul = target_Corg_total*(O2_U_clamp - target_O2_U_min)/(target_O2_U_max - target_O2_U_min)

        return Corg_cumul
    end

    # inverse of target_Corg_cumul
    function target_O2_U(Corg_cumul)
        Corg_cumul_clamp_norm = clamp(Corg_cumul/target_Corg_total, 0.0, 1.0)
        O2_U = target_O2_U_min + Corg_cumul_clamp_norm*(target_O2_U_max - target_O2_U_min)

        return O2_U
    end

    # calculate an approximation to the target Corg cumulative burial,
    # by adjusting burial in shelves by burial_adj_fac
    df_all_O2_U.burial_adj_fac = NaN*ones(nrow(df_all_O2_U))
    df_all_O2_U.Corg_b_cumsum_adj = NaN*ones(nrow(df_all_O2_U))

    local last_Corg_b_cumsum_adj = 0.0
    for irow in 1:nrow(df_all_O2_U)
        i = df_all_O2_U[irow, :i]
        ishelf = findfirst(x->x==i, ishelfX_floor) # nothing if i isn't a shelf floor box

        ocean_O2_U = df_all_O2_U[irow, :ocean_O2_U] # original version

        # attempt to centre bins on cumulative Corg
        # inextrow = min(irow + 1, nrow(df_all_O2_U))        
        # ocean_O2_U = 0.5*(df_all_O2_U[irow, :ocean_O2_U] + df_all_O2_U[inextrow, :ocean_O2_U]) 

        target_Corg_cl = target_Corg_cumul(ocean_O2_U)
       
        target_Corg_b = max(target_Corg_cl - last_Corg_b_cumsum_adj, 0.0)

        # attempt to centre bins on cumulative Corg
        # target_Corg_b = 2*target_Corg_b        

        # check in bounds
        target_Corg_b = min(target_Corg_b, target_Corg_adjust - last_Corg_b_cumsum_adj)
        target_Corg_b = max(target_Corg_b, 0.0)
   

        local Corg_b = df_all_O2_U[irow, :Corg_b]

        if isnothing(ishelf)
            # a "fixed" ocean box - leave it as is
            df_all_O2_U.Corg_b_cumsum_adj[irow] = last_Corg_b_cumsum_adj + Corg_b
        else
            # adjust to target
            adj_fac = target_Corg_b / Corg_b
            df_all_O2_U.burial_adj_fac[irow] = adj_fac
            df_all_O2_U.Corg_b_cumsum_adj[irow] = last_Corg_b_cumsum_adj + adj_fac*Corg_b
        end

        last_Corg_b_cumsum_adj = df_all_O2_U.Corg_b_cumsum_adj[irow]
    end

    # collect up shelf areas and adjustment factors
    shelf_Afloor, burial_adj_fac = Float64[], Float64[]
    for i in ishelfX_floor
        push!(shelf_Afloor, only(df_all_O2_U[df_all_O2_U.i .== i, :oceanfloor_Afloor]))
        push!(burial_adj_fac, only(df_all_O2_U[df_all_O2_U.i .== i, :burial_adj_fac]))
    end
    shelf_Afloor_adj = max.(shelf_Afloor.*burial_adj_fac, 1e6)
    P_O_romglb_shelf_Afloor_adj[fileroot] = (shelf_Afloor_adj, Corg_b_all, Corg_b_shelves)

    plot!(
        p_shelf_area_summary,
        shelf_Afloor_adj;
        label="target ($target_O2_U_min, $target_O2_U_max)",
        color=coloridx,
        seriestype=:stepmid,
    )

    p = plot(
        df_all_O2_U.Corg_b_cumsum, df_all_O2_U.ocean_O2_U;
        xlabel="Accumulated Corg burial (mol yr-1)", ylabel="frac O2 utilization", label="all",
        bottom_margin = 1Plots.mm, right_margin = 1Plots.mm,
    )
    plot!(
        p,
        df_all_O2_U.Corg_b_cumsum_adj, df_all_O2_U.ocean_O2_U;
        label="adjusted",
    )
    plot!(
        p,
        range(0.0, target_Corg_total, 100), target_O2_U;
        label="target",
    )
    display(p)

    plot!(
        p_corg_cumul_summary,
        df_all_O2_U.Corg_b_cumsum_adj, df_all_O2_U.ocean_O2_U;
        label="adj ($target_O2_U_min, $target_O2_U_max)",
        color=coloridx,
    )
    plot!(
        p_corg_cumul_summary,
        range(0.0, target_Corg_total, 100), target_O2_U;
        label="target ($target_O2_U_min, $target_O2_U_max)",
        color=coloridx,
        linestyle=:dash,
    )

    global coloridx += 1
end

println()
println("Adjusted shelf areas:")
for (fileroot, (v, corgb_tot, corgb_shelves)) in P_O_romglb_shelf_Afloor_adj
    println("fileroot: ", fileroot)
    println("unadjusted corg burial tot: ", corgb_tot, "     shelves: ", corgb_shelves)
    println("shelf_Afloor_adj: ", v)
    println()
end

plot!(
    p_shelf_area_summary;
    xlim=(0, 21), ylim=(2e11, 5e13),
    legend=:none,
)
display(p_shelf_area_summary)
savefig(p_shelf_area_summary, joinpath(output_figures_dir, "Shelf_Area_Summary.svg"))

plot!(
    p_corg_cumul_summary,
    ylim=(0.0, 1.0),
    legend=:none,
)
display(p_corg_cumul_summary)
savefig(p_corg_cumul_summary, joinpath(output_figures_dir, "Accumulated_Sum_Corg_burial.svg"))

# plot_vs_depth(paleorun.output, "fluxOceanBurial.flux_Corg",    "Corg burial (mol/yr)", "Depth (m)", NaN, "", true)
# @info "Global production=$(sum(bioprod[end])) (mol/yr), $(sum(bioprod[end])*12/1e15) (PgC/yr)"
