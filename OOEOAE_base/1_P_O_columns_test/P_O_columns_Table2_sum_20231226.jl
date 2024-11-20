using Logging
import DataFrames

import Plots

import PALEOboxes as PB
import PALEOmodel
# import PALEOreactions
import PALEOcopse

import Infiltrator
using Plots

include("../ReactionsOOEOAE_dev.jl")
# include("../SolverFunctionsOOEOAE.jl")
include("../SolverFunctionsOOEOAE2.jl")
include("../ooeoae_expts.jl")
# include("../ooeoae_plots.jl")

# dropbox_output_dir = "/Users/liziheng/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/1_P_O_columns_test"
# dropbox_output_dir = "C:/Users/ASUS/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/1_P_O_columns_test"
# dropbox_output_dir = "/home/sd336/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/1_P_O_columns_test"
dropbox_output_dir = joinpath(@__DIR__, "../../figures/1_P_O_columns_test")
output_figures_dir = joinpath(dropbox_output_dir, "P_O_columns_Table2_sum_20231226")

#####################################################
# Create model
#####################################################
model = PB.create_model_from_config(joinpath(@__DIR__, "P_O_columns.yaml"), "model1")

#####################################################
# Set experiment (parameters)
#####################################################
# basic setup
ooeoae_expts(
    model, [
        ("CPsea", 115.4, 461.6),
        ("land_flux_Bergman", false),
        ("k_anox", 100.0),  # slighly less sharp transition: default 1000 causes problems with fold detection
        ("anoxic_threshold", 0.5), # we want an estimate of corg-weighted-burial-flux under anoxic conditions
        # ("OPinit", (1.2*3.7e19, 2*3.1e15)),
        ("OPinit", (1.2*3.7e19, 2.5*3.1e15)),
    ]
)

# Experiments as listed in Table 2 in README
# All the tests in the Table 2 use preCambrian_Bergman's parameters
expts_table = [
    
    ####### 20231120, new plots for Fig2 ########
    ("PO 1", "modern_Bergman_lowO2U",        [("land_flux_Bergman", true), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 1.0)]),
    ("PO 2", "preCambrian_Bergman_lowO2U",        [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 1.0)]),
    ("PO 3", "preCambrian_Corgbf0p5_Bergman_lowO2U",        [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 0.5)]),

    ("PO 4", "modern_Bergman_neutrally_stable",     [("land_flux_Bergman", true), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.15, 0.6)), ("corg_burial_fac", 1.0)]),  
    ("PO 5", "preCambrian_Bergman_neutrally_stable",     [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.15, 0.6)), ("corg_burial_fac", 1.0)]),  
    ("PO 6", "preCambrian_Corgbf0p5_Bergman_neutrally_stable",     [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.15, 0.6)), ("corg_burial_fac", 0.5)]),  

    ("PO 7", "modern_Bergman_sharp_switch_unstable5",     [("land_flux_Bergman", true), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 1.0)]),  
    ("PO 8", "preCambrian_Bergman_sharp_switch_unstable5",     [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 1.0)]),  
    ("PO 9", "preCambrian_Corgbf0p5_Bergman_sharp_switch_unstable5",     [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 0.5)]),  

]


isdir(output_figures_dir) || mkdir(output_figures_dir)

P_O_columns_Table2_sum = Dict() # all results, indexed by fileroot

for (expt_id, fileroot, vector_pars) in expts_table
   # (fileroot, vector_pars) = expts_table[1]

   @info """

        ###############################################################################################################
        expt_id $expt_id   fileroot $fileroot
        ###############################################################################################################
     """

    ooeoae_expts(
        model, vector_pars
    )

    initial_state, modeldata = PALEOmodel.initialize!(model)

    # Find timeseries, P nullcline (critical manifold), O nullcline
    tspan = (0.0, 1e8)
    (paleorun, _, O_ts, P_ts) = SolverFunctionsOOEOAE2.find_time_series(model, initial_state, modeldata; has_A=false, tspan)

    pond = SolverFunctionsOOEOAE2.POnormDeriv(model, modeldata; fixed_values=initial_state, include_jacobian=true)
    # d/dP_m(dP_n/dt) of Jacobian = 0 at fold points
    pond_jac_P_P(P_n, O_n) = SolverFunctionsOOEOAE2.jacobian(pond, P_n, O_n)[1, 1]

    # P_lims = (0.0, 4.5)
    # O_lims = (0.0, 2.5)
    P_lims = (0.0, 5.5)
    O_lims = (0.0, 1.5)
    (;dPdt_line, dOdt_line, P_grid, O_grid) = SolverFunctionsOOEOAE2.find_nullclines_PO(
        pond; 
        P_lims, O_lims, P_npts=200, O_npts=200,
    )
    
    # eqb point where dO/dt = 0 as we move along dP/dt = 0 nullcline
    eqb_point = SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n)->pond(P_n, O_n)[2], dPdt_line; verbose=true)
    @info "eqb_point POA: $eqb_point"
    
    # fold point where d/dP_n(dP_n/dt) = 0 as we move along dP_n/dt=0 nullcline
    fold_points = SolverFunctionsOOEOAE2.Isoline.find_zeros_line(pond_jac_P_P, dPdt_line; verbose=true)
    @info "fold_points: $fold_points"

    # find the periodic
    t_ts = PB.get_data(paleorun.output, "global.tforce")

    # P burial and anoxia vs P at constant O
    mopb_const_O=Float64[]
    ANOX_const_O=Float64[]
    O_n = 0.5
    # @Infiltrator.infiltrate
    for P_n in P_grid
        (P_sms, O_sms) = pond(P_n, O_n)
        local burial_flux_total_P = only(PALEOmodel.get_array(pond.modeldata, "fluxOceanBurial.flux_total_P").values)
        push!(mopb_const_O, burial_flux_total_P)
        local AANOX = only(PALEOmodel.get_array(pond.modeldata, "oceanfloor.anoxia_burial_frac").values)
        push!(ANOX_const_O, AANOX) # approach 1
    end

    # anoxia over a grid
    ANOX_grid = fill(NaN, (length(O_grid), length(P_grid)))
    for (i, O_n) in enumerate(O_grid)
        for (j, P_n) in enumerate(P_grid)
            (P_sms, O_sms) = pond(P_n, O_n)
            local AANOX = only(PALEOmodel.get_array(pond.modeldata, "oceanfloor.anoxia_burial_frac").values)
            ANOX_grid[i, j] = AANOX
        end
    end
 
    # @Infiltrator.infiltrate

    results_namedtuple = (; 
        expt_id, paleorun, 
        t_ts, O_ts, P_ts, 
        dPdt_line, dOdt_line, P_grid, O_grid, eqb_point, fold_points,
        mopb_const_O, ANOX_const_O, ANOX_grid,
    )
 
    P_O_columns_Table2_sum[fileroot] = results_namedtuple

end

#################################
# summary figures 2
#################################
# | phase |_TS| phase |_TS|
# | stable|TS | sharp |TS |
# |_______|___|_______|___|
# |O2_U vs| local_| Pb vs.|
# |AccCorg| anoxia| P_norm|
# |_______|_______|_______|

linestyles=[:solid, :dash, :dashdotdot, :dashdot]

linewidths=[1, 1, 2, 2]

colors=[:green, :orange, :blue]

#### bottom left the O2_U vs. number of shelves (or accumulated Corg burial normalized)
#### bottom mid Corg weighted anoxia
#### bottom right the Pburial vs. P_norm
p_O2_U = Plots.plot(; 
    xlabel="Cumulative Corg burial", ylabel="O2 Utilization", xlims=(0, Inf), ylims=(0,1), legend=:none,
)
p_anoxia=Plots.plot(;
    xlabel="P_norm", ylabel="anoxia_frac", xlims=(0, Inf), ylims=(0, Inf), legend=:none,
)
ppp = Plots.plot(;
    ylabel="P_norm", xlabel="P burial (mol/yr)", ylims=(0, Inf), xlims=(0, Inf), legend=:none,
)

accumulated_Corg_burial_norm=collect(range(start=0, stop=1, length=100)) # number of shelves = 100

for (j, exptroot) in enumerate(["Bergman_lowO2U", "Bergman_neutrally_stable", "Bergman_sharp_switch_unstable5"])
    (; expt_id, paleorun, eqb_point, fold_points, P_grid, mopb_const_O, ANOX_const_O) =
        P_O_columns_Table2_sum["preCambrian_"*exptroot]
    O2_U = PB.get_data(paleorun.output, "ocean.O2_U")
    Plots.plot!(
        p_O2_U, accumulated_Corg_burial_norm, O2_U[end];
        color=:blue, linestyle=linestyles[j], linewidth=linewidths[j], label=expt_id, 
    )
    Plots.plot!(
        p_anoxia, P_grid, ANOX_const_O;
        color=:blue, linestyle=linestyles[j], linewidth=linewidths[j], label=expt_id,
    )
    Plots.plot!(
        ppp, mopb_const_O, P_grid;
        color=:blue, linestyle=linestyles[j], linewidth=linewidths[j], label=expt_id,       
    )
end


#### plot phase planes with 3x O nullclines #######
function plot_phase_plane(exptroot; linestyle_P_nullcline=:solid, linewidth_P_nullcline=1)
    p=Plots.plot(;
        xlims=[0.0, 5.0], ylims=[0.0, 1.5],          
        xlabel="P_norm", ylabel="O (PAL)", legend=:none,
    )

    # high CorgP x 3 nullclines
    # modern_Bergman preCambrian_Bergman preCambrian_Bergman_Corgbf0p5
    (;  expt_id, paleorun, 
        t_ts, O_ts, P_ts, dPdt_line, dOdt_line, 
        P_grid, O_grid, eqb_point, fold_points, mopb_const_O, ANOX_const_O, ANOX_grid,
    ) = P_O_columns_Table2_sum["preCambrian_$(exptroot)"]

    Pbal_P, Pbal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_line)
    # Obal_P, Obal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_line)
    Plots.plot!(
        p, Pbal_P, Pbal_O;
        color=:blue, linestyle=linestyle_P_nullcline, linewidth=linewidth_P_nullcline, label=false, 
    )
    (; element_counts, start_point_index, end_point_index) = SolverFunctionsOOEOAE2.find_periodic(
        P_ts, O_ts, t_ts, eqb_point[1];
        omit_t_start=5e6, # to estimate mid-point 
        include_t_start_crossings=false, # omit initial transient from zero crossings
        # plot=true,
    ) 
    if (end_point_index > 0) && (start_point_index > 0) # if the case is stable you may find the end_point_index==NaN
        Plots.plot!(
            p, P_ts[1:end_point_index], O_ts[1:end_point_index];
            color=:green, linestyle=:dash, label=false
        )
    else
        Plots.plot!(
            p, P_ts, O_ts;
            color=:green, linestyle=:dash, label=false
        )
    end

    Plots.contour!(
        p, P_grid, O_grid, ANOX_grid;
        levels=[1e-6, 0.25, 0.5, 0.75, 1.0-1e-6], color=:grey75,
        # levels=[1e-6, 0.5, 1.0], color=:black,
    )

    for (i, prefix) in enumerate(["modern_", "preCambrian_", "preCambrian_Corgbf0p5_"])
        local (; expt_id, paleorun, t_ts, O_ts, P_ts, dPdt_line, dOdt_line, P_grid, O_grid, eqb_point, fold_points, mopb_const_O, ANOX_const_O) =
            P_O_columns_Table2_sum[prefix*exptroot]
        # Pbal_P, Pbal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_line)
        Obal_P, Obal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_line)
        Plots.plot!(p, Obal_P, Obal_O, color=:red, linestyle=linestyles[i], label=expt_id)

        if !isempty(eqb_point)
            Plots.scatter!(p, [only(eqb_point)[1]], [only(eqb_point)[2]]; color=:red, markershape=:circle, markerstrokewidth=0, label=nothing);
        end
        Plots.scatter!(p,  map(x->x[1], fold_points), map(x->x[2], fold_points); color=:blue, markershape=:circle, markerstrokewidth=0, label=nothing);
    end

    return p
end


# plot phase plane with single O nullcline
function plot_phase_plane_single_O(exptroot; ylims=[0.0, 1.5], linestyle=:solid, linewidth=1, single_cycle=false)
    p_single_O=Plots.plot()
        
    (; expt_id, t_ts, O_ts, P_ts, dPdt_line, dOdt_line, eqb_point, fold_points) = P_O_columns_Table2_sum[exptroot]

    @info "exptroot = $exptroot, eqb_point = $eqb_point"
    (; element_counts, start_point_index, end_point_index) = SolverFunctionsOOEOAE2.find_periodic(
        P_ts, O_ts, t_ts, eqb_point[1]
    ) 
    println()

    Pbal_P, Pbal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_line)
    Obal_P, Obal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_line)
    p_single_O = Plots.plot(Pbal_P, Pbal_O; color=:blue, xlims=[0.0, 5.0], ylims, linestyle=linestyle, linewidth=linewidth, label=expt_id, xlabel="P_norm", ylabel="O (PAL)")
    Plots.plot!(p_single_O, Obal_P, Obal_O, color=:red, linestyle=:dash, label=false)

    if single_cycle == false
        Plots.plot!(p_single_O, P_ts, O_ts, color=:green, linestyle=linestyle, label=false)
    else
        Plots.plot!(p_single_O, P_ts[1:end_point_index], O_ts[1:end_point_index], color=:green, linestyle=linestyle, label=false)
    end

    Plots.scatter!(p_single_O, [only(eqb_point)[1]], [only(eqb_point)[2]]; color=:red, markershape=:circle, markerstrokewidth=0, label=nothing)
    Plots.scatter!(p_single_O,  map(x->x[1], fold_points), map(x->x[2], fold_points); color=:blue, markershape=:circle, markerstrokewidth=0, label=nothing);

    return p_single_O
end


##### time series for the waveform ######
function plot_time_series_waveform(
    exptroot;
    xlims=(0.0, 5.0e7),
    linestyle=:solid,
)
    p_time_series = Plots.plot(; xlims, xlabel="tmodel (yr)", ylims=(0, 4), ylabel="normalized value")

    (; t_ts, O_ts, P_ts, ) = P_O_columns_Table2_sum[exptroot]

    Plots.plot!(
        p_time_series, t_ts, P_ts; 
        color=:blue, linestyle, label="P_norm",
    )

    Plots.plot!(
        p_time_series, t_ts, O_ts; 
        color=:red, linestyle, label="O_norm",
    )

    return p_time_series
end


p_sum=Plots.plot()

# l = @layout[ a b c d;
#    e f g
# ]

l = @layout[
    a b c;
    d e f
]

p_sum = Plots.plot(
    plot_phase_plane("Bergman_lowO2U"), 
        plot_phase_plane("Bergman_neutrally_stable", linestyle_P_nullcline=:dash, linewidth_P_nullcline=1),
            plot_phase_plane("Bergman_sharp_switch_unstable5", linestyle_P_nullcline=:dashdotdot, linewidth_P_nullcline=2),
                # plot_phase_plane_single_O("preCambrian_Bergman_sharp_switch_unstable5"; linestyle=:dashdotdot, single_cycle=true, linewidth=2), 
               
    p_O2_U,
        # p_anoxia, 
        ppp,
            plot_time_series_waveform("preCambrian_Bergman_sharp_switch_unstable5"; linestyle=:dash),

    layout = l, 
    # size=(1200, 700),
    size=(1200, 650),
    left_margin = 5Plots.mm, right_margin = 5Plots.mm, bottom_margin = 5Plots.mm,
)

display(p_sum)

Plots.savefig(p_sum,  joinpath(output_figures_dir, "PO_secular_stability_ocillation_summary_20231226.svg"))
