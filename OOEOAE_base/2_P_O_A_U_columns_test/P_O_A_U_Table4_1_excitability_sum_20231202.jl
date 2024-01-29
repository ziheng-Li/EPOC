using Logging
import DataFrames
import Interpolations

using Plots
import GLMakie # import not using so doesn't conflict with Plots
# import MarchingCubes
# import GeometryBasics
# using Roots

import PALEOboxes as PB
import PALEOmodel
# import PALEOreactions
import PALEOcopse

# using Interpolations

# ReactionsOOEOAE_dev add Cisotopes, link DIC_sms -= local_Corgburial!
include("../ReactionsOOEOAE_dev.jl")
include("../SolverFunctionsOOEOAE2.jl")

include("../CarbBurial_dev.jl")
include("../Uranium.jl")

include("../ooeoae_expts.jl")
# include("../ooeoae_plots.jl")

# dropbox_output_dir = "/Users/liziheng/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/2_P_O_A_U_columns_test"
# dropbox_output_dir = "C:/Users/ASUS/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/2_P_O_A_U_columns_test"
# dropbox_output_dir = "C:/Users/sd336/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/1_P_O_columns_test"
dropbox_output_dir = joinpath(@__DIR__, "../../figures")
output_figures_dir = joinpath(dropbox_output_dir, "P_O_A_U_Table4_1_excitability_sum")

######################################
# key parameters that control plot appearance
include("expt_plot3D.jl")

#####################################################
# simplified dA/dt = oxidw + ocdeg + carbw + ccdeg - mocb - mccb
# -locb = 0, no -sfw
# modified after ../OOEOAE_ZHL/6_OOEOAE_columns_copsereload_Pw_test
#####################################################

model = PB.create_model_from_config(
    joinpath(@__DIR__, "../P_O_A_U_columns.yaml"), 
    "model1", 
    modelpars=Dict("CGconstant"=>false), # , "Aconstant"=>true
)

###############################################
# Set experiment (parameters)
###############################################
# basic setup 
ooeoae_expts(
    model, [
        ("k_O2_U", (0.23, 0.25)), 
        ("COPSE_locb", 1000.0, 0.0),
        ("CPsea", 115.4, 461.6), # base setup values for Bergman_COPSE
        ("k_anox", 100.0),
        # ("OPAinit", (0.5*3.7e19, 4*3.1e15, 3*3.193e18)),
        ("corg_burial_fac", 1.25), # Corgb2 = 1.2 too close to the folder, changed to 1.35
    ]
)

# Experiments as listed in Table2 in README
expts_table = [
    # reproduce the test case (same name) in Table2, as the blank run
    # ("preCambrian_Bergman_unstable_basline1_Psilw_only",            [("P_weathering", 3.9e10, 1.0, 0.0, 0.0),   ("COPSE_locb", 1000.0, 0.0),    ("k_O2_U", (0.30, 0.70))]), # marginally unstable
    
    # ############# anoxic + decrease degassing #############
    # ## excitability test, same duration, two magnitude of pulse
    ##########################################################

    # # the duration of the CO2pulse [P_cycle_time, A_cycle_time] 
    # this is a OOE case
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb7_Psilw_only_CO2pulse1", [("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
    #      ("CO2pulse", [5e7, 5e7+1.0e6], [-2.0e12, -2.0e12]), ("OPAinit", (0.591*3.7e19, 3.75*3.1e15, 3.94*3.193e18))]), 
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb7_Psilw_only_CO2pulse2", [("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
    #      ("CO2pulse", [5e7, 5e7+1.0e6], [-3.0e12, -3.0e12]), ("OPAinit", (0.591*3.7e19, 3.75*3.1e15, 3.94*3.193e18))]), 

    # the duration of the CO2pulse < P_cycle_time 
    # this is a OAE case
    ("POA 14", "preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_CO2pulse1", [("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
         ("CO2pulse", [0.1e7, 0.1e7+1.0e4], [1.0e14, 1.0e14]), ("OPAinit", (0.31*3.7e19, 1.0*3.1e15, 4.07*3.193e18))]),
    ("POA 15", "preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_CO2pulse2", [("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
         ("CO2pulse", [0.1e7, 0.1e7+1.0e4], [1.6e14, 1.6e14]), ("OPAinit", (0.31*3.7e19, 1.0*3.1e15, 4.07*3.193e18))]),

]

isdir(output_figures_dir) || mkdir(output_figures_dir)

P_O_A_U_columns_Table4_1 = Dict()

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

    ################################
    # get time_series
    ####################################
    (paleorun, A_ts, O_ts, P_ts) = SolverFunctionsOOEOAE2.find_time_series(model, initial_state, modeldata; has_A=true, tspan)

    ###########################################################
    # phase plane plot with nullclines and folds
    ########################################################
    include_jacobian = true
    # poand: initial_state, poand_end: end_state
    poand = SolverFunctionsOOEOAE2.POAnormDeriv(model, modeldata; fixed_values=initial_state, fixed_t=tspan[1], include_jacobian)
    poand_end = SolverFunctionsOOEOAE2.POAnormDeriv(model, modeldata; fixed_values=initial_state, fixed_t=tspan[end], include_jacobian)
    if include_jacobian
        # d/dP_m(dP_n/dt) component of Jacobian to identify folds
        # -ve is an attracting surface, +ve is repelling
        global poand_jac_P_P(P_n, O_n, A_n) = SolverFunctionsOOEOAE2.jacobian(poand, P_n, O_n, A_n)[1, 1]
        global poand_jac_P_P_end(P_n, O_n, A_n) = SolverFunctionsOOEOAE2.jacobian(poand_end, P_n, O_n, A_n)[1, 1]
    else
        # label everything as 'attracting'
        global poand_jac_P_P(P_n, O_n, A_n) = -1.0 
        global poand_jac_P_P_end(P_n, O_n, A_n) = -1.0 
    end

    # get nullcline surfaces
    (;dPdt_surf, dOdt_surf, dAdt_surf, A_grid) = SolverFunctionsOOEOAE2.find_nullclines_POA(poand; P_lims, O_lims, A_lims=(2.0, A_lims[2]))
    (dPdt_surf_end, dOdt_surf_end, dAdt_surf_end, _, _, A_grid_end) = SolverFunctionsOOEOAE2.find_nullclines_POA(poand_end; P_lims, O_lims, A_lims=(2.0, A_lims[2]))
   
 
    ########################################################
    # dA/dt = 0 surface, for initial_state
    dAdt_surf_P, dAdt_surf_O, dAdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dAdt_surf)
    # dA/dt = 0 intersection with dP/dt = 0   
    dAdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dPdt_surf))
    # split at folds
    dAdt_dPdt_linesegs_POA = SolverFunctionsOOEOAE2.split_line_sign_f(dAdt_dPdt_line_POA, poand_jac_P_P)

    # do the same thing as above, but for the end_state!
    dAdt_surf_P_end, dAdt_surf_O_end, dAdt_surf_A_end = SolverFunctionsOOEOAE2.tuples_to_coords(dAdt_surf_end)
    dAdt_dPdt_line_POA_end = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand_end(P_n, O_n, A_n)[3], dPdt_surf_end))
    dAdt_dPdt_linesegs_POA_end = SolverFunctionsOOEOAE2.split_line_sign_f(dAdt_dPdt_line_POA_end, poand_jac_P_P_end)
    ########################################################

    ########################################################
    # dO/dt = 0 surface
    dOdt_surf_P, dOdt_surf_O, dOdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_surf)
    # dO/dt = 0 intersection with dP/dt = 0
    dOdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[2], dPdt_surf))
    # split at folds
    dOdt_dPdt_linesegs_POA = SolverFunctionsOOEOAE2.split_line_sign_f(dOdt_dPdt_line_POA, poand_jac_P_P)
    ########################################################
       
    # find folds (P component of Jacobian = 0)
    folds_dPdt_lines_POA = SolverFunctionsOOEOAE2.Isoline.find_isolines(poand_jac_P_P, dPdt_surf)
    # order folds so first index is lowest P
    if length(folds_dPdt_lines_POA) == 2
        f_1, f_2 = folds_dPdt_lines_POA[1], folds_dPdt_lines_POA[2]
        P_1, P_2 = first(f_1)[1], first(f_2)[1]
        if P_1 > P_2
            folds_dPdt_lines_POA[1:2] .= f_2, f_1
        end
    end

    A_target_val = 4.0 # roughly midpoint of limit cycle
    A_ks = findfirst(x -> x > A_target_val, A_grid)
    A_val = A_grid[A_ks] # actual A used
    dPdt_linesegs_const_A = SolverFunctionsOOEOAE2.split_line_sign_f(dPdt_surf[:, A_ks], poand_jac_P_P)
    Obal_P, Obal_O, _ = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_surf[:, 1]) # dO/dt = 0 is independent of A so just pick first A

    O_val = 0.32 # roughly midpoint of limit cycle
    dPdt_line_const_O, dAdt_line_const_O = SolverFunctionsOOEOAE2.find_nullclines_const_O(poand, O_val; P_lims, A_lims)
    # split at folds
    dPdt_linesegs_const_O = SolverFunctionsOOEOAE2.split_line_sign_f(dPdt_line_const_O, poand_jac_P_P)

    # find where dA/dt=0 along line with dO/dt=0 and dP/dt=0
    eqb_point = only(SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dOdt_dPdt_line_POA; verbose=true))
    @info "eqb_point POA: $eqb_point"

    ##########################################
    # save everything in 3D for the sum plot!
    P_O_A_U_columns_Table4_1[fileroot] = 
        (; expt_id, paleorun, modeldata, t_ts, A_ts, P_ts, O_ts, A_val, O_val,
            dPdt_surf, A_grid, dAdt_surf, dOdt_surf,
            dPdt_linesegs_const_A, Obal_P, Obal_O, 
            dPdt_line_const_O, dAdt_line_const_O, dPdt_linesegs_const_O,
            folds_dPdt_lines_POA, dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA_end, eqb_point)
    ##########################################

end # end of for loop


#####################

# linestyles = [:solid, :dash, :dashdot, :dot, :dashdotdot]

####### 3D phase plane #######
for (j, exptroot) in enumerate(["CO2pulse1", "CO2pulse2"])
    local (; expt_id, paleorun, modeldata, t_ts, A_ts, P_ts, O_ts, A_val, O_val,
        dPdt_surf, A_grid, dAdt_surf, dOdt_surf,
        dPdt_linesegs_const_A, Obal_P, Obal_O, 
        dPdt_line_const_O, dAdt_line_const_O, dPdt_linesegs_const_O,
        folds_dPdt_lines_POA, dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA_end, eqb_point) =
        P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_"*exptroot]
    
    fig, axs, pltobj = plot_3D_critical_manifold(dPdt_surf, A_grid, folds_dPdt_lines_POA)

    # add intersection lines of O_nullcline and slow-manifold surface
    SolverFunctionsOOEOAE2.plot_segments!(
        GLMakie.lines!, axs, dOdt_dPdt_linesegs_POA, [1, 2, 3], (:dash, :solid), nothing; 
        color=:red, linewidth=glm_linewidth,
    )
    # add intersection lines of A_nullcline and slow-manifold surface
    SolverFunctionsOOEOAE2.plot_segments!(
        GLMakie.lines!, axs, dAdt_dPdt_linesegs_POA, [1, 2, 3], (:dash, :solid), nothing;
        color=:black, linewidth=glm_linewidth,
    )
    # add limit cycle time series
    (; element_counts, start_point_index, end_point_index) = SolverFunctionsOOEOAE2.find_periodic(P_ts, O_ts, t_ts, eqb_point) 
    GLMakie.lines!(
        axs,
        P_ts[start_point_index:end_point_index], O_ts[start_point_index:end_point_index], A_ts[start_point_index:end_point_index];
        color=:green, linestyle=:solid, linewidth=2*glm_linewidth,
    )

    GLMakie.scatter!(axs, eqb_point; color=:red, markersize=glm_eqb_markersize,)

    GLMakie.save(joinpath(output_figures_dir, "preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_"*exptroot*"_3D.png"), fig;
    px_per_unit=5, # Makie 0.20 increase resolution of saved figure (600 x 5 = 3000 pixels)
    ) # closes window
end

# ####### 2D plots in single one ########
(; expt_id, paleorun, modeldata, t_ts, A_ts, P_ts, O_ts, A_val, O_val,
        dPdt_surf, A_grid, dAdt_surf, dOdt_surf,
        dPdt_linesegs_const_A, Obal_P, Obal_O, 
        dPdt_line_const_O, dAdt_line_const_O, dPdt_linesegs_const_O,
        folds_dPdt_lines_POA, dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA_end, eqb_point)  = 
        P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_CO2pulse1"]

plot_PA = SolverFunctionsOOEOAE2.plot_3D_mapping_to_PA_phase(dPdt_linesegs_const_O, dAdt_line_const_O, Obal_P, Obal_O, O_val, P_ts, A_ts, true, true, [0.0, 4.0], [2.0, 5.0])
plot_OA = SolverFunctionsOOEOAE2.plot_3D_mapping_to_OA_phase(folds_dPdt_lines_POA, dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA_end, O_ts, A_ts, true, true, [0.0, 1.0], [2.0, 5.0], linestyles[1])

for (i, exptroot) in enumerate(["CO2pulse2"])
    local (; expt_id, paleorun, modeldata, t_ts, A_ts, P_ts, O_ts, A_val, O_val,
    dPdt_surf, A_grid, dAdt_surf, dOdt_surf,
    dPdt_linesegs_const_A, Obal_P, Obal_O, 
    dPdt_line_const_O, dAdt_line_const_O, dPdt_linesegs_const_O,
    folds_dPdt_lines_POA, dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA_end, eqb_point)  = 
        P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_"*exptroot]
    Plots.plot!(plot_OA, O_ts, A_ts; color=:green, linestyles=linestyles[i+1], label=false)
    Plots.scatter!(plot_OA, [eqb_point[2]], [eqb_point[3]]; color=:red, markershape=:circle, markerstrokewidth=0, label=nothing);

    Plots.plot!(plot_PA, P_ts, A_ts; color=:green, linestyles=linestyles[i+1], label=false)
    Plots.scatter!(plot_PA, [eqb_point[1]], [eqb_point[3]]; color=:red, markershape=:circle, markerstrokewidth=0, label=nothing);
end

############# time series ###############
plot_d13C = plot()
for (i, exptroot) in enumerate(["CO2pulse1", "CO2pulse2"])
    (; paleorun) = 
    P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_"*exptroot]
    ts = PB.get_data(paleorun.output, "global.tforce")
    DIC_delta = PB.get_data(paleorun.output, "ocean.DIC_delta")
    Plots.plot!(plot_d13C, ts, DIC_delta; color=:green, linestyles=linestyles[i], xlabel="tmodel (yr)", ylabel="d13C (‰)", label=false, xlims=(0e7,4e6))
end

plot_ANOX = plot()
for (i, exptroot) in enumerate(["CO2pulse1", "CO2pulse2"])
    (; paleorun) = 
        P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_"*exptroot]
    ts = PB.get_data(paleorun.output, "global.tforce")
    anoxia_burial_frac = PB.get_data(paleorun.output, "oceanfloor.anoxia_burial_frac")
    Plots.plot!(plot_ANOX, ts, anoxia_burial_frac; color=:green, linestyles=linestyles[i], xlabel="tmodel (yr)", ylabel="ANOX", label=false, xlims=(0e7,4e6))
end

plot_2_1 = plot()
for (i, exptroot) in enumerate(["CO2pulse1", "CO2pulse2"])
    (; paleorun, P_ts) = 
    P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_"*exptroot]
    ts = PB.get_data(paleorun.output, "global.tforce")
    Plots.plot!(plot_2_1, ts, P_ts; color=:green, linestyles=linestyles[i], ylabel="P_norm", label=exptroot, xlims=(0e7,4e6))
end

plot_2_2 = plot()
for (i, exptroot) in enumerate(["CO2pulse1", "CO2pulse2"])
    (; paleorun) = 
    P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_"*exptroot]
    ts = PB.get_data(paleorun.output, "global.tforce")
    CO2pulse_delta=PB.get_data(paleorun.output, "global.CO2pulse")
    CO2pulse=[]
    for i = 1:length(CO2pulse_delta)
        push!(CO2pulse, CO2pulse_delta[i].v)
    end
    Plots.plot!(plot_2_2, ts, CO2pulse; color=:green, linestyles=linestyles[i], ylabel="CO2pulse", label=false, xlims=(0e7,4e6))
end

l = @layout[
            grid(2,1){0.5w} grid(4,1)
]

p_sum = Plots.plot(plot_OA, plot_PA,
                    plot_2_2, plot_2_1, plot_d13C, plot_ANOX, layout = l, left_margin = 5Plots.mm, bottom_margin = 1Plots.mm, size=(800, 600))

savefig(p_sum, joinpath(output_figures_dir, "Excitability_CO2pulse_summary_20231202.svg"))

# # # report back for the initial OPA state
# # # for (i, exptroot) in enumerate(["Corgb4", "Corgb3", "Corgb6"])
# #     (_, _, A_ts, P_ts, O_ts, _, _, _, _, _,Obal_P, Obal_O, _, _, _, _, _, _) = 
# #     P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_CO2pulse1"]
# #     println(" A=$(A_ts[end]), P=$(P_ts[end]), O=$(O_ts[end])")
# # # end

