using Logging
import DataFrames

using Plots 
# using Roots

import PALEOboxes as PB
import PALEOmodel
# import PALEOreactions
import PALEOcopse

# ReactionsOOEOAE_dev add Cisotopes, link DIC_sms -= local_Corgburial!
include("../../../src/ReactionsOOEOAE_dev.jl")
include("../../../src/SolverFunctionsOOEOAE2.jl")

include("CarbBurial_dev.jl")
include("Uranium.jl")

include("../ooeoae_expts.jl")
# include("../ooeoae_plots.jl")

# Archived figures
# dropbox_output_dir = joinpath(@__DIR__, "../../figures/2_P_O_A_U_columns_test")

# Local figures
isdir("figures") || mkdir("figures")
dropbox_output_dir = "figures"

output_figures_dir = joinpath(dropbox_output_dir, "P_O_A_U_Table3_1")
isdir(output_figures_dir) || mkdir(output_figures_dir)

#####################################################
# simplified dA/dt = oxidw + ocdeg + carbw + ccdeg - mocb - mccb
# -locb = 0, no -sfw
# modified after ../OOEOAE_ZHL/6_OOEOAE_columns_copsereload_Pw_test
#####################################################

model = PB.create_model_from_config(
    joinpath(@__DIR__, "P_O_A_U_columns.yaml"), 
    "model1", 
    # modelpars=Dict("CGconstant"=>false), # , "Aconstant"=>true
    modelpars=Dict(
        "CGconstant"=>true, # fixed sedimentary reservoirs (small effect on degassing reductant input)
        "Aconstant"=>true, # fix A (should make no difference as no weathering fluxes used, prescribed P input flux)
    )
)

###############################################
# Set experiment (parameters)
###############################################
# basic setup 
# this is a P-O-A model configuration, but with land-surface P weathering removed to create
# what is effectively a P-O only model
# NB: there are two prescribed P inputs: P_pulse (used for time dependent weathering) and P_weathering (to set a background constant)

P_norm_val = 3.1e15 # mol (as set in yaml file)
O_norm_val = 3.7e19 # mol (as set in yaml file)
ooeoae_expts(
    model, [
        ("k_O2_U", (0.23, 0.25)), # sharp switch
        ("corg_burial_fac", 2.1), # stable oxic
        ("CPsea", 115.4, 461.6), # base setup values for Bergman_COPSE
        ("OPinit", (0.5*O_norm_val, 0.15*P_norm_val)),
        ("P_weathering", 39e9-20e9),  # background constant P weathering input (NB: land weathering P input already has default of zero)
    ]
)


# Excitability and rate-dependent ramping tests

# names for plot labelling
pulse_labels = Dict("Ppulse1"=>"ΔP_norm 0.24", "Ppulse2"=>"ΔP_norm 0.48") # total P input for excitability tests / normalized P

# rate_labels = Dict("rate1"=>"20 My", "rate2"=>"10 My", "rate3"=>"5 My", "rate4"=>"0.1 My")
rate_times = [20, 5, 1, 0.1]
rate_labels = Dict("rate1"=>"20 My", "rate2"=>"5 My", "rate3"=>"1 My", "rate4"=>"0.1 My")

# NB: this is a P-O-A model yaml config, but with land-weathering P input set to zero to create a P-O only model
# and then Ppulse is added to P_weathering
expts_table = [ 
    # excitation of sharp-switch case by a near-instantaneous P input pulse
    # # fileroot                                      k10_phosw  sil carb ox                    CPland   k11_landfrac
    ("PO 10", "modern_Bergman_sharp_switch_stable_oxic_Ppulse1", [("Ppulse", [2e7, 2e7+5e4], [1.5e10, 1.5e10*5e4/5e4])]),  # small pike
    ("PO 11", "modern_Bergman_sharp_switch_stable_oxic_Ppulse2", [("Ppulse", [2e7, 2e7+5e4], [3.0e10, 3.0e10*5e4/5e4])]),  # full round

    # ("modern_Bergman_sharp_switch_stable_oxic_Ppulse_rate0", [# ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
    #     ("Ppulse", [0, 1.5e7, 1.5e7+2.0e7, 1e8], [0, 0, 0, 0])]),  # no Pw feedback

    # ramp of P input, for globally stable case
    ("PO 12", "modern_Bergman_lowO2U_Ppulse_rate1", [# ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
        ("Ppulse", [0, 1.5e7, 1.5e7+1e6*rate_times[1], 1e8], [0, 0, 20e9, 20e9]), ("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 2.1)]),  # no Pw feedback

    ("PO 13", "modern_Bergman_lowO2U_Ppulse_rate2", [# ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
        ("Ppulse", [0, 1.5e7, 1.5e7+1e6*rate_times[2], 1e8], [0, 0, 20e9, 20e9]), ("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 2.1)]),  # no Pw feedback
    
    ("PO 14", "modern_Bergman_lowO2U_Ppulse_rate3", [# ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
        ("Ppulse", [0, 1.5e7, 1.5e7+1e6*rate_times[3], 1e8], [0, 0, 20e9, 20e9]), ("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 2.1)]),  # no Pw feedback
    
    ("PO 15", "modern_Bergman_lowO2U_Ppulse_rate4", [# ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
        ("Ppulse", [0, 1.5e7, 1.5e7+1e6*rate_times[4], 1e8], [0, 0, 20e9, 20e9]), ("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 2.1)]),  # no Pw feedback

    # ramp of P input, for sharp switch case
    ("PO 16", "modern_Bergman_sharp_switch_stable_oxic_Ppulse_rate1", [# ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
        ("Ppulse", [0, 1.5e7, 1.5e7+1e6*rate_times[1], 1e8], [0, 0, 20e9, 20e9]), ("k_O2_U", (0.23, 0.25)),]),  # no Pw feedback

    ("PO 17", "modern_Bergman_sharp_switch_stable_oxic_Ppulse_rate2", [# ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
        ("Ppulse", [0, 1.5e7, 1.5e7+1e6*rate_times[2], 1e8], [0, 0, 20e9, 20e9]), ("k_O2_U", (0.23, 0.25)),]),  # no Pw feedback
    
    ("PO 18", "modern_Bergman_sharp_switch_stable_oxic_Ppulse_rate3", [# ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
        ("Ppulse", [0, 1.5e7, 1.5e7+1e6*rate_times[3], 1e8], [0, 0, 20e9, 20e9]), ("k_O2_U", (0.23, 0.25)),]),  # no Pw feedback
    
    ("PO 129", "modern_Bergman_sharp_switch_stable_oxic_Ppulse_rate4", [# ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
        ("Ppulse", [0, 1.5e7, 1.5e7+1e6*rate_times[4], 1e8], [0, 0, 20e9, 20e9]), ("k_O2_U", (0.23, 0.25)),]),  # no Pw feedback

]


P_O_A_U_columns_table3_1 = Dict() # all results, indexed by fileroot

for (expt_id, fileroot, vector_pars) in expts_table
    # (fileroot, vector_pars) = expts_table[1]

    ooeoae_expts(
        model, vector_pars
    )

    tspan = (0.0, 1e8) # yr # tspan=(-1000e6, 0)

    #########################################################
    # Initialize
    #########################################################

    initial_state, modeldata = PALEOmodel.initialize!(model)

    (paleorun, _, O_ts, P_ts) = SolverFunctionsOOEOAE2.find_time_series(model, initial_state, modeldata; has_A=false, tspan)

    pond = SolverFunctionsOOEOAE2.POnormDeriv(model, modeldata; fixed_values=initial_state, include_jacobian=true)
    pond_end = SolverFunctionsOOEOAE2.POnormDeriv(model, modeldata; fixed_values=initial_state, fixed_t=tspan[end], include_jacobian=true)
    # d/dP_m(dP_n/dt) of Jacobian = 0 at fold points
    pond_jac_P_P(P_n, O_n) = SolverFunctionsOOEOAE2.jacobian(pond, P_n, O_n)[1, 1]

    P_lims = (0.0, 4.5)
    O_lims = (0.0, 2.5)
    (;dPdt_line, dOdt_line, P_grid, O_grid) = SolverFunctionsOOEOAE2.find_nullclines_PO(pond; P_lims, O_lims)
    (dPdt_line_end, dOdt_line_end, P_grid_end, O_grid_end) = SolverFunctionsOOEOAE2.find_nullclines_PO(pond_end; P_lims, O_lims)

    # eqb point where dO/dt = 0 as we move along dP/dt = 0 nullcline
    eqb_point = SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n)->pond(P_n, O_n)[2], dPdt_line; verbose=true)
    @info "eqb_point: $eqb_point"
    eqb_point_end = SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n)->pond_end(P_n, O_n)[2], dPdt_line_end; verbose=true)
        @info "eqb_point_end POA: $eqb_point_end"

    # fold point where d/dP_n(dP_n/dt) = 0 as we move along dP_n/dt=0 nullcline
    fold_points = SolverFunctionsOOEOAE2.Isoline.find_zeros_line(pond_jac_P_P, dPdt_line; verbose=true)
    @info "fold_points: $fold_points"

    # find the periodic
    t_ts = PB.get_data(paleorun.output, "global.tforce")
    # (periodic, start_point_index, end_point_index) = SolverFunctionsOOEOAE2.find_periodic(P_ts, O_ts, t_ts) 

    results_namedtuple = (; expt_id, paleorun, t_ts, O_ts, P_ts, dPdt_line, dOdt_line, P_grid, O_grid, dPdt_line_end, dOdt_line_end, P_grid_end, O_grid_end, eqb_point, eqb_point_end, fold_points)

    P_O_A_U_columns_table3_1[fileroot] = results_namedtuple
    Pbal_P, Pbal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_line)
    Obal_P, Obal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_line)

    # ######### Plot #########
    # gr(size=(1000, 600)) # plotlyjs(size=(1000, 600))
    # pager=PALEOmodel.PlotPager(
    #     (3, 2), (legend_background_color=nothing, );
    #     displayfunc=(plot, nplot)->savefig(plot, joinpath(output_figures_dir, "$(fileroot).svg")), # save to file instead of default display
    # )
    # pager(
    #     (
    #         Plots.plot(Pbal_P, Pbal_O, xlims=[0.0, 4.5], ylims=[0.0, 2.5], label="slow manifold", xlabel="P_norm", ylabel="O (PAL)");
    #         Plots.plot!(Pbal_P, Pbal_O, xlims=[0.0, 4.5], ylims=[0.0, 2.5], label="slow manifold", xlabel="P_norm", ylabel="O (PAL)");
    #         Plots.plot!(Obal_P, Obal_O, label="O nullcline");
    #         Plots.plot!(P_ts, O_ts, label="time series", left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
    #     ),

    #     Plots.plot(ylabel="Reservoirs", paleorun.output, ["atmocean.O_norm", "ocean.P_norm", "sedcrust.C_norm", "sedcrust.G_norm", "atmocean.A_norm", "ocean.U_norm"]),
    #     # plot(paleorun.output, ylabel="P fluxes", ["fluxRtoOcean.flux_P", "fluxOceanBurial.flux_total_P"]), 
    #     plot(paleorun.output, ylabel="O fluxes", ["fluxOceanBurial.flux_total_Corg", "fluxAtoLand.flux_O2", "fluxSedCrusttoAOcean.flux_Redox"]), 

    #     Plots.plot(ylabel="TEMP (K)", paleorun.output, ["global.TEMP"]),
    #     Plots.plot(title="Carbon isotopes",  paleorun.output, ["atmocean.A_delta", "ocean.DIC_delta", "atm.CO2_delta"], ylabel="delta 13C (per mil)"), # , xlims=(-1e6, 10e6), "ocean.mccb_delta", "sedcrust.C_delta" ; extrakwargs...),
    #     Plots.plot(title="Uranium isotopes",  paleorun.output, ["fluxRtoOcean.flux_U.v_delta", "ocean.U.v_delta"], ylabel="d238U/235U"),
    #     # plot(title="anoxia_burial_frac", paleorun.output, ["oceanfloor.anoxia_burial_frac"])
    # )

    # pager(:newpage)
    
    # plot_anim_P_O(Pbal_P,Pbal_O,Obal_P,Obal_O,P_ts,O_ts,"plot_table3_1/$(fileroot)_anim.gif")
end

# #################################
# # summary figures (single panels)
# #################################
linestyles=[:solid, :dash, :dashdot, :dashdotdot]
linewidths=[2, 1, 1, 1]

function plot_phaseplane_multi_trajectory(
    root::String,
    post_fix::Vector{String};
    plot_P_nullcline_end=true,  # false to omit 2nd P nullcline
    plot_fold_points=true, # false to omit fold points
    trajectory_labels=Dict(l=>l for l in post_fix),
    ts_start=15e6, # first time to use (eg to omit spinup)
    xlims=(0,5),
    ylims=(0,2.25),
    label_P_nullcline="P nullcline initial",
    label_P_nullcline_end="P nullcline final",
    label_O_nullcline="O nullcline"
)
    p2 = plot(; 
        xlims=xlims, ylims=ylims,
        xlabel="P_norm", ylabel="O (PAL)",
    )

    (; t_ts, O_ts, P_ts, dPdt_line, dOdt_line, dPdt_line_end, eqb_point, eqb_point_end, fold_points) =
        # P_O_A_U_columns_table3_1["modern_Bergman_lowO2U_Ppulse_rate1"] 
        P_O_A_U_columns_table3_1[root*post_fix[1]]

    Pbal_P, Pbal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_line)   
    plot!(p2, Pbal_P, Pbal_O; color=:blue, label=label_P_nullcline)

    if plot_P_nullcline_end
        Pbal_P_end, Pbal_O_end = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_line_end)
        plot!(p2, Pbal_P_end, Pbal_O_end; color=:blue, linestyle=:dash, label=label_P_nullcline_end)
    end

    Obal_P, Obal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_line)
    plot!(p2, Obal_P, Obal_O, color=:red; label=label_O_nullcline)

    Plots.scatter!(p2, [only(eqb_point)[1]], [only(eqb_point)[2]]; color=:red, markershape=:circle, markerstrokewidth=0, label=nothing)
    Plots.scatter!(p2, [only(eqb_point_end)[1]], [only(eqb_point_end)[2]]; color=:red, markershape=:circle, markerstrokewidth=0, label=nothing)
    if plot_fold_points
        Plots.scatter!(p2,  map(x->x[1], fold_points), map(x->x[2], fold_points); color=:blue, markershape=:circle, markerstrokewidth=0, label=nothing)
    end

    for (exptroot, ls, lw)  in zip(post_fix, linestyles, linewidths)
        (; t_ts, O_ts, P_ts) =
            P_O_A_U_columns_table3_1[root*exptroot]
        # omit times < ts_start
        i_first = findfirst(x->x >= ts_start, t_ts)
        plot!(
            p2, P_ts[i_first:end], O_ts[i_first:end];
            label=trajectory_labels[exptroot], color=:green, linestyle=ls, linewidth=lw,
        )
    end

    return p2
end

function plot_multi_timeseries_P_norm(
    root::String,
    post_fix::Vector{String};
    trajectory_labels=Dict(l=>l for l in post_fix),
    xlims=(0e7, 3.5e7),
    # ylims=(0,2.25),
)
    p1_1 = plot(
        xlabel="tmodel (yr)",  xlims=xlims, 
        ylabel="P_norm",
        right_margin = 5Plots.mm,
    )
    for (expt, ls)  in zip(post_fix, linestyles)
        (; t_ts, P_ts) =
            P_O_A_U_columns_table3_1[root*expt] 

        plot!(
            p1_1, t_ts, P_ts;
            color=:green, linestyle=ls, label=trajectory_labels[expt],
        )
    end

    return p1_1
end

function plot_multi_timeseries_P_input(
    root::String,
    post_fix::Vector{String};
    trajectory_labels=Dict(l=>l for l in post_fix),
    xlims=(0e7, 3.5e7),
    # ylims=(0,2.25),
)
    p1_1 = plot(
        xlabel="tmodel (yr)",  xlims=xlims, 
        ylabel="P input (1e10 mol yr-1)",
        right_margin = 5Plots.mm,
    )
    for (expt, ls)  in zip(post_fix, linestyles)
        (; paleorun, t_ts) =
            P_O_A_U_columns_table3_1[root*expt] 
      
       P_input = (
                PB.get_data(paleorun.output, "fluxRtoOcean.flux_P") # background constant P_weathering
            .+  PB.get_data(paleorun.output, "global.Ppulse") # Ppulse is applied directly to ocean.P_sms
        )
     
        plot!(
            p1_1, t_ts, P_input./1e10;
            color=:green, linestyle=ls, label=trajectory_labels[expt],
        )
    end

    return p1_1
end


# plot options
ts_xlims = (10e6, 35e6) # time series xlims (P perturbation starts at 15e6 yr)


l = @layout[
    a b c
]
p_sum = Plots.plot(
    plot_phaseplane_multi_trajectory(
        "modern_Bergman_sharp_switch_stable_oxic_", ["Ppulse1", "Ppulse2"];
        trajectory_labels=pulse_labels,
        xlims=(0, 2.1), ylims=(0, 0.75),
        plot_P_nullcline_end=false,
        label_P_nullcline=nothing, label_P_nullcline_end=nothing, label_O_nullcline=nothing,
    ),
        plot_phaseplane_multi_trajectory(
            "modern_Bergman_lowO2U_Ppulse_", ["rate1", "rate2", "rate3", "rate4"];
            trajectory_labels=rate_labels,
            xlims=(0, 4.2), ylims=(0, 1.5),
            plot_fold_points=false,
            label_P_nullcline=nothing, label_P_nullcline_end=nothing, label_O_nullcline=nothing,
        ),
            plot_phaseplane_multi_trajectory(
                "modern_Bergman_sharp_switch_stable_oxic_Ppulse_", ["rate1", "rate2", "rate3", "rate4"];
                trajectory_labels=rate_labels,
                xlims=(0, 4.2), ylims=(0, 1.5),
                label_P_nullcline=nothing, label_P_nullcline_end=nothing, label_O_nullcline=nothing,
            ),
            
    layout = l, 
    left_margin = 5Plots.mm, right_margin = 5Plots.mm, bottom_margin = 5Plots.mm,
    size=(900, 300),
)

display(p_sum)
savefig(p_sum, joinpath(output_figures_dir, "P_O_excitability_rate_dependent_oxic_20231227.svg"))


l = @layout[
    grid(3,3)
]

p_sum_SI = Plots.plot(
    plot_phaseplane_multi_trajectory(
        "modern_Bergman_sharp_switch_stable_oxic_", ["Ppulse1"];
        trajectory_labels=pulse_labels,
        xlims=(0, 2.1), ylims=(0, 0.75),
        plot_P_nullcline_end=false,
        label_P_nullcline="P nullcline",
    ),
        plot_phaseplane_multi_trajectory(
            "modern_Bergman_lowO2U_Ppulse_", ["rate1", "rate2", "rate3", "rate4"];
            trajectory_labels=rate_labels,
            xlims=(0, 4.2), ylims=(0, 1.5),
            plot_fold_points=false,
        ),
            plot_phaseplane_multi_trajectory(
                "modern_Bergman_sharp_switch_stable_oxic_Ppulse_", ["rate1", "rate2", "rate3", "rate4"];
                trajectory_labels=rate_labels,
                xlims=(0, 4.2), ylims=(0, 1.5)
            ),
            

    plot_phaseplane_multi_trajectory(
        "modern_Bergman_sharp_switch_stable_oxic_", ["Ppulse2"];
        trajectory_labels=pulse_labels,
        xlims=(0, 2.1), ylims=(0, 0.75),
        plot_P_nullcline_end=false,
        label_P_nullcline="P nullcline",
    ),
        plot_multi_timeseries_P_input(
            "modern_Bergman_lowO2U_Ppulse_", ["rate1", "rate2", "rate3", "rate4"];
            trajectory_labels=rate_labels,
            xlims=ts_xlims,
        ),
            plot_multi_timeseries_P_input(
                "modern_Bergman_sharp_switch_stable_oxic_Ppulse_", ["rate1", "rate2", "rate3", "rate4"];
                trajectory_labels=rate_labels,
                xlims=ts_xlims,
            ),
           

    plot_multi_timeseries_P_norm(
        "modern_Bergman_sharp_switch_stable_oxic_", ["Ppulse1", "Ppulse2"];
        trajectory_labels=pulse_labels,
        xlims=ts_xlims,
    ),
        plot_multi_timeseries_P_norm(
            "modern_Bergman_lowO2U_Ppulse_", ["rate1", "rate2", "rate3", "rate4"];
            trajectory_labels=rate_labels,
            xlims=ts_xlims,
        ),
            plot_multi_timeseries_P_norm(
                "modern_Bergman_sharp_switch_stable_oxic_Ppulse_", ["rate1", "rate2", "rate3", "rate4"];
                trajectory_labels=rate_labels,
                xlims=ts_xlims,
            ),
            
    layout = l, 
    left_margin = 5Plots.mm, right_margin = 5Plots.mm, bottom_margin = 1Plots.mm,
    size=(1200, 1000),
)

display(p_sum_SI)
savefig(p_sum_SI, joinpath(output_figures_dir, "P_O_excitability_rate_dependent_oxic_SI_20231227.svg"))


