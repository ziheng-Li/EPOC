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
include("../../../src/ReactionsOOEOAE_dev.jl")
include("../../../src/SolverFunctionsOOEOAE2.jl")

include("CarbBurial_dev.jl")
include("Uranium.jl")

include("../ooeoae_expts.jl")


# Archived figures
# dropbox_output_dir = joinpath(@__DIR__, "../../figures/2_P_O_A_U_columns_test")

# Local figures
isdir("figures") || mkdir("figures")
dropbox_output_dir = "figures"

output_figures_dir = joinpath(dropbox_output_dir, "P_O_A_U_Table2_3D_sum_graphicstest")
isdir(output_figures_dir) || mkdir(output_figures_dir)

######################################
# key parameters that control plot appearance
include("expt_plot3D.jl")

##########################################

#####################################################
# simplified dA/dt = oxidw + ocdeg + carbw + ccdeg - mocb - mccb
# -locb = 0, no -sfw
# modified after ../OOEOAE_ZHL/6_OOEOAE_columns_copsereload_Pw_test
#####################################################

model = PB.create_model_from_config(
    joinpath(@__DIR__, "P_O_A_U_columns.yaml"), 
    "model1", 
    modelpars=Dict("CGconstant"=>false), # , "Aconstant"=>true
)

###############################################
# Set experiment (parameters)
###############################################
# basic setup 
ooeoae_expts(
    model, [
        # ("land_flux_Bergman", false),
        ("k_O2_U", (0.23, 0.25)), # sharp switch
        ("CPsea", 115.4, 461.6), # base setup values for Bergman_COPSE
        ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
        ("CO2pulse", [0, 5e7, 7.0e7, 1e8], [0, 0, 0, 0]), 
        ("COPSE_locb", 1000.0, 0.0),
        ("k_anox", 100.0),
    ]
)

# Experiments as listed in Table2 in README
expts_table = [
    # test case for the corg_burial_fac (Onulcline)
    # # stable case
    # ("preCambrian_Bergman_Psilw_only", [("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 0.8), ("OPAinit", (1.0*3.7e19, 0.0*3.1e15, 1*3.193e18))]),
    ("preCambrian_Bergman_Psilw_only", [("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 0.8), ("OPAinit", (0.1*3.7e19, 2.0*3.1e15, 5.0*3.193e18))]),

    # corg_burial_fac=[1.131, 0.476]
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb1_Psilw_only", [("corg_burial_fac", 1.5)]),
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only", [("corg_burial_fac", 1.2)]),
    ("preCambrian_Bergman_sharp_switch_unstable5_Corgb3_Psilw_only", [("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 1.1), ("OPAinit", (0.32*3.7e19, 1.06*3.1e15, 4.19*3.193e18))]),
    ("preCambrian_Bergman_sharp_switch_unstable5_Corgb4_Psilw_only", [("k_O2_U", (0.23, 0.25)),("corg_burial_fac", 0.7), ("OPAinit", (0.548*3.7e19, 1.21*3.1e15, 4.53*3.193e18))]),
    # # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb5_Psilw_only", [("corg_burial_fac", 0.6)]),
    ("preCambrian_Bergman_sharp_switch_unstable5_Corgb6_Psilw_only", [("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 0.5), ("OPAinit", (0.59*3.7e19, 3.53*3.1e15, 3.83*3.193e18))]),
    # # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb7_Psilw_only", [("corg_burial_fac", 0.47)]),
    # # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb8_Psilw_only", [("corg_burial_fac", 0.4)]),

    # # corg_burial_fac=[1.111, 0.708]
    # # sharpest case with low CPsea, lossing the PO dominated region
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb1_lowCPsea_Psilw_only", [("CPsea", 115.4, 230.8), ("corg_burial_fac", 1.5)]),
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb2_lowCPsea_Psilw_only", [("CPsea", 115.4, 230.8), ("corg_burial_fac", 1.2)]),
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb3_lowCPsea_Psilw_only", [("CPsea", 115.4, 230.8), ("corg_burial_fac", 1.0)]),
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb4_lowCPsea_Psilw_only", [("CPsea", 115.4, 230.8), ("corg_burial_fac", 0.8)]),
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb5_lowCPsea_Psilw_only", [("CPsea", 115.4, 230.8), ("corg_burial_fac", 0.6)]),
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb6_lowCPsea_Psilw_only", [("CPsea", 115.4, 230.8), ("corg_burial_fac", 0.5)]),
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb7_lowCPsea_Psilw_only", [("CPsea", 115.4, 230.8), ("corg_burial_fac", 0.47)]),
    # ("preCambrian_Bergman_sharp_switch_unstable5_Corgb8_lowCPsea_Psilw_only", [("CPsea", 115.4, 230.8), ("corg_burial_fac", 0.4)]),

    # # corg_burial_fac=[1.066, 0.371]
    # ("preCambrian_Corgbf0p5_Bergman_unstable_baseline2_lowO2U_Corgb1_Psilw_only",  [("k_O2_U", (0.125, 0.25)), ("corg_burial_fac", 1.5)]),
    # ("preCambrian_Corgbf0p5_Bergman_unstable_baseline2_lowO2U_Corgb2_Psilw_only",  [("k_O2_U", (0.125, 0.25)), ("corg_burial_fac", 1.2)]),
    # ("preCambrian_Corgbf0p5_Bergman_unstable_baseline2_lowO2U_Corgb3_Psilw_only",  [("k_O2_U", (0.125, 0.25)), ("corg_burial_fac", 1.0)]),
    # ("preCambrian_Corgbf0p5_Bergman_unstable_baseline2_lowO2U_Corgb4_Psilw_only",  [("k_O2_U", (0.125, 0.25)), ("corg_burial_fac", 0.8)]),
    # ("preCambrian_Corgbf0p5_Bergman_unstable_baseline2_lowO2U_Corgb5_Psilw_only",  [("k_O2_U", (0.125, 0.25)), ("corg_burial_fac", 0.6)]),
    # ("preCambrian_Corgbf0p5_Bergman_unstable_baseline2_lowO2U_Corgb6_Psilw_only",  [("k_O2_U", (0.125, 0.25)), ("corg_burial_fac", 0.5)]),
    # ("preCambrian_Corgbf0p5_Bergman_unstable_baseline2_lowO2U_Corgb7_Psilw_only",  [("k_O2_U", (0.125, 0.25)), ("corg_burial_fac", 0.4)]),
    # ("preCambrian_Corgbf0p5_Bergman_unstable_baseline2_lowO2U_Corgb8_Psilw_only",  [("k_O2_U", (0.125, 0.25)), ("corg_burial_fac", 0.3)]),

]

P_O_A_U_columns_Table2_3D_sum = Dict()

for (fileroot, vector_pars) in expts_table

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
    t_ts = PB.get_data(paleorun.output, "global.tforce")

    ###########################################################
    # phase plane plot with nullclines and folds
    ########################################################
    include_jacobian = true
    poand = SolverFunctionsOOEOAE2.POAnormDeriv(model, modeldata; include_jacobian)
    if include_jacobian
        # d/dP_m(dP_n/dt) component of Jacobian to identify folds
        # -ve is an attracting surface, +ve is repelling
        global poand_jac_P_P(P_n, O_n, A_n) = SolverFunctionsOOEOAE2.jacobian(poand, P_n, O_n, A_n)[1, 1]
    else
        # label everything as 'attracting'
        global poand_jac_P_P(P_n, O_n, A_n) = -1.0 
    end

    
    # get nullcline surfaces
    (;dPdt_surf, dOdt_surf, dAdt_surf, A_grid) = SolverFunctionsOOEOAE2.find_nullclines_POA(
        poand; P_lims, O_lims, A_lims=(2.0, A_lims[2])
    )
   

    # dPdt surface
    dPdt_surf_P, dPdt_surf_O, dPdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_surf)
 
    # dA/dt = 0 surface
    dAdt_surf_P, dAdt_surf_O, dAdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dAdt_surf)
    # dA/dt = 0 intersection with dP/dt = 0   
    dAdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dPdt_surf))
    # split at folds
    dAdt_dPdt_linesegs_POA = SolverFunctionsOOEOAE2.split_line_sign_f(dAdt_dPdt_line_POA, poand_jac_P_P)
 
    # dO/dt = 0 surface
    dOdt_surf_P, dOdt_surf_O, dOdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_surf)
    # dO/dt = 0 intersection with dP/dt = 0
    dOdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[2], dPdt_surf))
    # split at folds
    dOdt_dPdt_linesegs_POA = SolverFunctionsOOEOAE2.split_line_sign_f(dOdt_dPdt_line_POA, poand_jac_P_P)
       
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

    ##########################################
    # 2D mappings
    ############################################

    Plots.gr(size=(1400, 800)) # plotlyjs(size=(1600, 900))
    pager=PALEOmodel.PlotPager(
        (2, 3),
        (legend_background_color=nothing,);
        displayfunc=(plot, nplot)->savefig(plot, joinpath(output_figures_dir, "$(fileroot).svg")),  # save to file instead of default display)
    )

    # P, O phase plane at constant A: find dP/dt = 0 for a representative value of A
    # (can get this from surfaces we have already calculated)
    A_target_val = 4.0 # roughly midpoint of limit cycle
    A_ks = findfirst(x -> x > A_target_val, A_grid)
    A_val = A_grid[A_ks] # actual A used
    dPdt_linesegs_const_A = SolverFunctionsOOEOAE2.split_line_sign_f(dPdt_surf[:, A_ks], poand_jac_P_P)
    Obal_P, Obal_O, _ = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_surf[:, 1]) # dO/dt = 0 is independent of A so just pick first A
    # plot_PO = SolverFunctionsOOEOAE2.plot_3D_mapping_to_PO_phase(dPdt_linesegs_const_A, Obal_P, Obal_O, A_val, P_ts, O_ts, [0.0, 5.0], [0.0, 1.5])

    # P, A plane at constant O: find dP/dt = 0, dA/dt = 0 for one value of O
    # (need to calculate nullclines in P, A surface at constant O)
    O_val = 1.5 # roughly midpoint of limit cycle
    dPdt_line_const_O, dAdt_line_const_O = SolverFunctionsOOEOAE2.find_nullclines_const_O(poand, O_val; P_lims, A_lims)
    # split at folds
    dPdt_linesegs_const_O = SolverFunctionsOOEOAE2.split_line_sign_f(dPdt_line_const_O, poand_jac_P_P)
    # plot_PA = SolverFunctionsOOEOAE2.plot_3D_mapping_to_PA_phase(dPdt_linesegs_const_O, dAdt_line_const_O, O_val, P_ts, A_ts, [0.0, 5.0], [2.0, 5.0])

    # O, A plane, "projecting" critical manifold ie using P value(s) from critical manifold at each O, A
    dAdt_dPdt_linesegs_POA_end = Vector{Any}()
    # plot_OA = SolverFunctionsOOEOAE2.plot_3D_mapping_to_OA_phase(folds_dPdt_lines_POA, dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA_end, O_ts, A_ts, [0.0, 1.5], [2.0, 5.0])

    # find where dA/dt=0 along line with dO/dt=0 and dP/dt=0
    eqb_point = only(SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dOdt_dPdt_line_POA; verbose=true))
    @info "eqb_point POA: $eqb_point"

    # find the periodic
    ((periodic, start_point_index, end_point_index)) = SolverFunctionsOOEOAE2.find_periodic(P_ts, O_ts, t_ts, eqb_point) 

    ##########################################
    # save everything in 3D for the sum plot!
    P_O_A_U_columns_Table2_3D_sum[fileroot] = 
        (; paleorun, modeldata, t_ts, A_ts, P_ts, O_ts, A_val, O_val,
            dPdt_surf, A_grid, dAdt_surf, dOdt_surf,
            dPdt_linesegs_const_A, Obal_P, Obal_O, 
            dPdt_line_const_O, dAdt_line_const_O, dPdt_linesegs_const_O,
            folds_dPdt_lines_POA, dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA_end, eqb_point)
    ##########################################

end # end of for loop



# ####### Figure3 of the draft ########
# # | Stable| 3D   |
# # |  3D   |sharp |
# # |_______|______|
# # |_TS_|_TS_|_TS_|

# # ####### For the 3D part #######

# linestyles = [:solid, :dash, :dashdot, :dot, :dashdotdot]
linestyles = [:solid, :solid, :solid]
# linewidths=[1, 1, 2, 2]

# One panel with "central" limit cycle
(; dPdt_surf, A_grid, folds_dPdt_lines_POA, t_ts, A_ts, P_ts, O_ts,
    dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, eqb_point) = 
        P_O_A_U_columns_Table2_3D_sum["preCambrian_Bergman_sharp_switch_unstable5_Corgb4_Psilw_only"]

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
(periodic, start_point_index, end_point_index) = SolverFunctionsOOEOAE2.find_periodic(P_ts, O_ts, t_ts, eqb_point) 
GLMakie.lines!(
    axs,
    P_ts[start_point_index:end_point_index], O_ts[start_point_index:end_point_index], A_ts[start_point_index:end_point_index];
    color=:green, linestyle=:solid, linewidth=2*glm_linewidth,
)

GLMakie.scatter!(axs, eqb_point; color=:red, markersize=glm_eqb_markersize,)
display(fig)
GLMakie.save(joinpath(output_figures_dir, "POAU_figure3_sharp_switch_3D_phase_plane_central.png"), fig; px_per_unit=5,)
# redisplay needed to reset scaling ?
display(fig)

# One panel with two "marginally unstable" time series
# dP/dt=0 is the same so use any one of 3
(; dPdt_surf, A_grid, folds_dPdt_lines_POA,) = 
        P_O_A_U_columns_Table2_3D_sum["preCambrian_Bergman_sharp_switch_unstable5_Corgb4_Psilw_only"]

fig, axs, pltobj = plot_3D_critical_manifold(dPdt_surf, A_grid, folds_dPdt_lines_POA)

for (i, exptroot) in zip([2, 3], ["Corgb3", "Corgb6"]) # ["Corgb4", "Corgb3", "Corgb6"]) # 
    local (; paleorun, t_ts, A_ts, P_ts, O_ts,
        dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, eqb_point) = 
        P_O_A_U_columns_Table2_3D_sum["preCambrian_Bergman_sharp_switch_unstable5_"*exptroot*"_Psilw_only"]

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
    (periodic, start_point_index, end_point_index) = SolverFunctionsOOEOAE2.find_periodic(P_ts, O_ts, t_ts, eqb_point) 
    GLMakie.lines!(
        axs,
        P_ts[start_point_index:end_point_index], O_ts[start_point_index:end_point_index], A_ts[start_point_index:end_point_index];
        color=:green, linestyle=linestyles[i], linewidth=2*glm_linewidth,
    )

    GLMakie.scatter!(axs, eqb_point; color=:red, markersize=glm_eqb_markersize,)
end
display(fig)
GLMakie.save(joinpath(output_figures_dir, "POAU_figure3_sharp_switch_3D_phase_plane_marginal.png"), fig; px_per_unit=5,)
# redisplay needed to reset scaling ?
display(fig)

# # ####### For the 2D TS part #######
plot_time_series_PO = Dict()
for (i, exptroot) in enumerate(["Corgb4", "Corgb3", "Corgb6"]) # 
    local (; t_ts, A_ts, P_ts, O_ts) = 
        P_O_A_U_columns_Table2_3D_sum["preCambrian_Bergman_sharp_switch_unstable5_"*exptroot*"_Psilw_only"]

    
    # for time series
    pp = Plots.plot()
    Plots.plot!(pp, t_ts, O_ts, label="O_norm", linestyle=linestyles[i], color=:red)
    Plots.plot!(pp, t_ts, P_ts, label="P_norm", linestyle=linestyles[i], color=:blue)
    Plots.plot!(pp, t_ts, A_ts, label="A_norm", linestyle=linestyles[i], color=:black)
    plot_time_series_PO[exptroot] = pp
end
l = @layout[a b c]
p_sum = Plots.plot(plot_time_series_PO["Corgb3"], plot_time_series_PO["Corgb4"], plot_time_series_PO["Corgb6"],
                     layout = l, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, size=(1200, 200))
savefig(p_sum, joinpath(output_figures_dir, "POAU_figure3_time_series.svg"))


####### Figure3 of the draft the stable case, show all the surface ########
(;paleorun, modeldata, t_ts, A_ts, P_ts, O_ts, A_val, O_val,
                dPdt_surf, dAdt_surf, dOdt_surf,
                dPdt_linesegs_const_A, Obal_P, Obal_O, 
                dPdt_line_const_O, dAdt_line_const_O, dPdt_linesegs_const_O,
                folds_dPdt_lines_POA, dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA_end, eqb_point) = 
        P_O_A_U_columns_Table2_3D_sum["preCambrian_Bergman_Psilw_only"]

fig, axs, pltobj = plot_3D_critical_manifold(dPdt_surf, A_grid, folds_dPdt_lines_POA)

# dA/dt = 0 surface
plot_3D_A_nullcline(axs, dAdt_surf)
SolverFunctionsOOEOAE2.plot_segments!(
    GLMakie.lines!, axs, dAdt_dPdt_linesegs_POA, [1, 2, 3], (:dash, :solid), nothing;
    color=:black, linewidth=glm_linewidth,
)

# dO/dt = 0 surface
plot_3D_O_nullcline(axs, dOdt_surf)
SolverFunctionsOOEOAE2.plot_segments!(
    GLMakie.lines!, axs, dOdt_dPdt_linesegs_POA, [1, 2, 3], (:dash, :solid), nothing;
    color=:red, linewidth=glm_linewidth,
)

# TIME seeries
GLMakie.lines!(axs, P_ts, O_ts, A_ts; color=:green, linewidth=2*glm_linewidth,)

GLMakie.scatter!(axs, eqb_point; color=:red, markersize=glm_eqb_markersize,)

display(fig)

GLMakie.save(joinpath(output_figures_dir, "POAU_figure3_stable_3D_phase_plane.png"), fig; px_per_unit=5,)
# redisplay needed to reset scaling ?
display(fig)


