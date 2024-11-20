using Logging
import DataFrames
import Plots

using Plots 
using Roots

import PALEOboxes as PB
import PALEOmodel
# import PALEOreactions
import PALEOcopse

include("../ReactionsOOEOAE_dev.jl")
# include("../SolverFunctionsOOEOAE.jl")
include("../SolverFunctionsOOEOAE2.jl")
include("../ooeoae_expts.jl")


# dropbox_output_dir = "/Users/liziheng/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/1_P_O_columns_test"
# dropbox_output_dir = "C:/Users/ASUS/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/1_P_O_columns_test"
# dropbox_output_dir = "/home/sd336/Dropbox/BACE_OOEOAE/DainesLiOverleaf/OOEOAE_base/1_P_O_columns_test"
dropbox_output_dir = joinpath(@__DIR__, "../../figures/1_P_O_columns_test")
output_figures_dir = joinpath(dropbox_output_dir, "P_O_columns_FigS1_20231202")


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
    ]
)

expts_table = [
    # fileroot, expt_id, params_vec
    # all these are experiment ID PO_XX 
    # modern cases with locb (high Corg burial)
    ("101",  "modern_Bergman_constCPsea1",         [("land_flux_Bergman", true), ("CPsea", 115.4, 115.4), ("k_O2_U", (0.15, 0.85)), ("corg_burial_fac", 1.0)]),
    ("104",  "modern_Bergman_constCPsea2",         [("land_flux_Bergman", true), ("CPsea", 461.6, 461.6), ("k_O2_U", (0.15, 0.85)), ("corg_burial_fac", 1.0)]),
    ("107",  "modern_Bergman_lowO2U",              [("land_flux_Bergman", true), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 1.0)]),
    ("1010", "modern_Bergman_lowCPsea",            [("land_flux_Bergman", true), ("CPsea", 115.4, 230.8), ("k_O2_U", (0.15, 0.85)), ("corg_burial_fac", 1.0)]),
    
    # preCambrian cases no locb (1.0 mocb)
    ("102", "preCambrian_Bergman_constCPsea1",     [("land_flux_Bergman", false), ("CPsea", 115.4, 115.4), ("k_O2_U", (0.15, 0.85)), ("corg_burial_fac", 1.0)]),
    ("105", "preCambrian_Bergman_constCPsea2",     [("land_flux_Bergman", false), ("CPsea", 461.6, 461.6), ("k_O2_U", (0.15, 0.85)), ("corg_burial_fac", 1.0)]),
    ("108",  "preCambrian_Bergman_lowO2U",         [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 1.0)]),
    ("1011", "preCambrian_Bergman_lowCPsea",       [("land_flux_Bergman", false), ("CPsea", 115.4, 230.8), ("k_O2_U", (0.15, 0.85)), ("corg_burial_fac", 1.0)]),
    
    # preCambrian cases no locb (0.5 mocb)
    ("103", "preCambrian_Corgbf0p5_Bergman_constCPsea1",    [("land_flux_Bergman", false), ("CPsea", 115.4, 115.4), ("k_O2_U", (0.15, 0.85)), ("corg_burial_fac", 0.5)]),
    ("106", "preCambrian_Corgbf0p5_Bergman_constCPsea2",     [("land_flux_Bergman", false), ("CPsea", 461.6, 461.6), ("k_O2_U", (0.15, 0.85)), ("corg_burial_fac", 0.5)]),
    ("109", "preCambrian_Corgbf0p5_Bergman_lowO2U",         [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.05, 0.75)), ("corg_burial_fac", 0.5)]),
    ("1012", "preCambrian_Corgbf0p5_Bergman_lowCPsea",        [("land_flux_Bergman", false), ("CPsea", 115.4, 230.8), ("k_O2_U", (0.15, 0.85)), ("corg_burial_fac", 0.5)]),
    ################################

    # fix O nullcline and eqb point at that of S26, adjust "sharpness"
    # eqb_point: [(2.412486059397197, 0.5613435384538765)]
    # Calculation here is:
    # "oxic" fold at P=1.0, (O - P*k_O2_U[2]) = 0.0 -> O = k_O2_U[2]
    # "anoxic" fold at P=4.0, (O - P*k_O2_U[1]) = 0.0 -> O = 4*k_O2_U[1]
    # so between folds, 
    #   O(P) = k_O2_U[2] + (4*k_O2_U[1] - k_O2_U[2])*(P - 1)/3
    # we require this line to include P_eqb, O_eqb, so
    #   0.561 = k_O2_U[2] + (4*k_O2_U[1] - k_O2_U[2])*(2.412 - 1)/3
    #   0.561 = k_O2_U[2] + 4*0.471*k_O2_U[1] - 0.471 * k_O2_U[2]
    #   0.561 = 4*0.471*k_O2_U[1] + (1 - 0.471)*k_O2_U[2]
    #   k_O2_U[1] = (0.561 - 0.529*k_O2_U[2])/(4*0.471)
    ("1013",  "preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness0", [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.132, 0.59)), ("corg_burial_fac", 0.65)]),
    ("1014",  "preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness4", [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.561/4, 0.561)), ("corg_burial_fac", 0.65)]),
    ("1015",  "preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness3", [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.145, 0.53)), ("corg_burial_fac", 0.65)]),
    ("1016", "preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness2",  [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.157, 0.5)), ("corg_burial_fac", 0.65)]), 
    ("1017", "preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness1",   [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.185, 0.4)), ("corg_burial_fac", 0.65)]), 
    ("1018",  "preCambrian_Bergman_sharp_switch_unstable5_waveform",           [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 0.65)]), 

    #################################
    # Bergman_sharp_switch_unstable5 wave form
    # share the case S18
    ("1019",  "modern_Bergman_sharp_switch_unstable5_waveform",                [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 1.0)]), 
    ("1020",  "preCambrian_Corgbf0p5_Bergman_sharp_switch_unstable5_waveform",  [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 0.48)]), 

    #################################
    # tested here not included in the paper!!
    # ("T1", "modern_Bergman_marginally_stable1",  [("land_flux_Bergman", true), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.1, 0.75)), ("corg_burial_fac", 1.0)]),
    # ("T2", "preCambrian_Bergman_marginally_stable1", [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.1, 0.75)), ("corg_burial_fac", 0.8)]),
    # ("T3", "preCambrian_Corgbf0p5_Bergman_marginally_stable1",  [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.1, 0.75)), ("corg_burial_fac", 0.4)]),
    # ("T4", "modern_Bergman_marginally_stable2",  [("land_flux_Bergman", true), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.18, 0.75)), ("corg_burial_fac", 1.0)]),
    # ("T5",  "preCambrian_Bergman_marginally_stable2",  [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.18, 0.75)), ("corg_burial_fac", 0.8)]),
    # ("T6", "preCambrian_Corgbf0p5_Bergman_marginally_stable2", [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.18, 0.75)), ("corg_burial_fac", 0.4)]),
    # ("T7", "modern_Bergman_unstable_baseline2",  [("land_flux_Bergman", true), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.25, 0.75)), ("corg_burial_fac", 1.0)]), 
    # ("T8", "preCambrian_Bergman_unstable_baseline2",   [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.25, 0.75)), ("corg_burial_fac", 0.8)]), 
    # ("T9", "preCambrian_Corgbf0p5_Bergman_unstable_baseline2", [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.25, 0.75)), ("corg_burial_fac", 0.4)]),
    # ("T10", "modern_Bergman_sharp_switch_unstable5", [("land_flux_Bergman", true), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 1.0)]), 
    # ("T11", "preCambrian_Bergman_sharp_switch_unstable5",   [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 0.8)]), 
    # ("T12", "preCambrian_Corgbf0p5_Bergman_sharp_switch_unstable5",  [("land_flux_Bergman", false), ("CPsea", 115.4, 461.6), ("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 0.4)]), 

    # #################################
    # # the following tests are for low CPsea sharpness cases 
    # # get the value here, used by the 3D phaseplane cases
    # # the low CPsea case in P_O system has corg_burial_fac=1.0, 
    # # the low CPsea case in P_O_A system has corg_burial_fac=0.8, 
    # # i choose the value 1.o here which is more closer to the mid of unstable regime !
    # # eqb_point: [(1.5497874316538034, 0.5455092964561753)]
    # # Calculation here is:
    # # "oxic" fold at P=1.0, (O - P*k_O2_U[2]) = 0.0 -> O = k_O2_U[2]
    # # "anoxic" fold at P=2.0, (O - P*k_O2_U[1]) = 0.0 -> O = 2*k_O2_U[1]
    # # # eqb_point: [(1.4180738546345013, 0.3361302216677178)]
    # # # -> k_O2_U2 = k_O2_U1*2 - (k_O2_U1*2 - 0.33613)/(2-1.41807)
    # # # -> k_O2_U1 = (k_O2_U2 - (0.33613)/(2-1.41807))/(2-2/(2-1.41807))
    # # ("T13", "preCambrian_Bergman_lowCPsea",       [("land_flux_Bergman", false), ("CPsea", 115.4, 230.8), ("k_O2_U", (0.15, 0.85)), ("corg_burial_fac", 0.9)]), # intersection point is not [(1.4180738546345013, 0.3361302216677178)]
    # ("T14", "preCambrian_Bergman_lowCPsea_sharpness1",       [("land_flux_Bergman", false), ("CPsea", 115.4, 230.8), ("k_O2_U", (0.123, 0.4)), ("corg_burial_fac", 0.9)]),
    # ("T15", "preCambrian_Bergman_lowCPsea_sharpness2",       [("land_flux_Bergman", false), ("CPsea", 115.4, 230.8), ("k_O2_U", (0.158, 0.35)), ("corg_burial_fac", 0.9)]),
    # ("T16", "preCambrian_Bergman_lowCPsea_sharpness3",       [("land_flux_Bergman", false), ("CPsea", 115.4, 230.8), ("k_O2_U", (0.336/2, 0.336)), ("corg_burial_fac", 0.9)]),
    # ("T17", "preCambrian_Bergman_lowCPsea_sharpness4",       [("land_flux_Bergman", false), ("CPsea", 115.4, 230.8), ("k_O2_U", (0.179, 0.32)), ("corg_burial_fac", 0.9)]),
    # ("T18", "preCambrian_Bergman_lowCPsea_sharpness5",       [("land_flux_Bergman", false), ("CPsea", 115.4, 230.8), ("k_O2_U", (0.193, 0.30)), ("corg_burial_fac", 0.9)]),
    # ("T19", "preCambrian_Bergman_lowCPsea_sharpness6",       [("land_flux_Bergman", false), ("CPsea", 115.4, 230.8), ("k_O2_U", (0.23, 0.25)), ("corg_burial_fac", 0.9)]),
]

isdir(output_figures_dir) || mkdir(output_figures_dir)

P_O_columns_FigS1 = Dict() # all results, indexed by fileroot


for (expt_id, fileroot, vector_pars) in expts_table
    # (expt_id, fileroot, vector_pars) = expts_table[1]

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
 
     P_lims = (0.0, 4.5)
     O_lims = (0.0, 8.0)
     (;dPdt_line, dOdt_line, P_grid, O_grid) = SolverFunctionsOOEOAE2.find_nullclines_PO(pond; P_lims, O_lims)
     
     # eqb point where dO/dt = 0 as we move along dP/dt = 0 nullcline
     eqb_point = SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n)->pond(P_n, O_n)[2], dPdt_line; verbose=true)
     @info "eqb_point: $eqb_point"

    # fold point where d/dP_n(dP_n/dt) = 0 as we move along dP_n/dt=0 nullcline
    fold_points = SolverFunctionsOOEOAE2.Isoline.find_zeros_line(pond_jac_P_P, dPdt_line; verbose=true)
    @info "fold_points: $fold_points"

    # find the periodic
    t_ts = PB.get_data(paleorun.output, "global.tforce")
    # (; element_counts, start_point_index, end_point_index, sign_change) = SolverFunctionsOOEOAE2.find_periodic(P_ts, O_ts, t_ts; Spec_P=P_ts[end])
 
    P_x = collect(range(0.0, 8.0, length=100))
    mopb=[]
    ANOX=[]
    O_n = 0.5
    P_w = PB.get_data(paleorun.output, "land.P_weath")
    for P_n in P_x
        (P_sms, O_sms) = pond(P_n, O_n)
        push!(mopb, P_w[end] - P_sms*3.1e15) # P_initial = 3.1e15
        local AANOX = PALEOmodel.get_array(pond.modeldata, "oceanfloor.anoxia_burial_frac").values
        push!(ANOX, AANOX[1]) # approach 1
    end
 
    results_namedtuple = (; expt_id, paleorun, t_ts, O_ts, P_ts, dPdt_line, dOdt_line, P_grid, O_grid, eqb_point, fold_points, mopb, ANOX)
 
    P_O_columns_FigS1[fileroot] = results_namedtuple
 
end

########## SI figure, support the main text figure 2
# |__|__|__|__| # phase plane and time series
# |__|__|__|__|
# |     |     | # left O2_U vs. acc_Corg
# |_____|_____| # right Pburial vs. P_norm

# plots_vector=Vector{Any}()

linestyles=[:solid, :dash, :dashdot, :dashdotdot, :solid, :dash, :dashdot, :dashdotdot, :solid, :dash, :dashdot, :dashdotdot, :solid, :dash, :dashdot, :dashdotdot]

colors=[:black, :black, :red, :red, :red, :red, :blue, :blue, :blue, :blue, :blue, :blue]

linewidths=[2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1]


# #### bottom left the O2_U vs. number of shelves (or accumulated Corg burial normalized)
# pp = plot()

# accumulated_Corg_burial_norm=collect(range(start=0, stop=1, length=100)) # number of shelves = 100

# for (j, exptroot) in enumerate(["Bergman_constCPsea1", "Bergman_constCPsea2", "Bergman_lowO2U", "Bergman_lowCPsea", "Bergman_marginally_stable1", "Bergman_marginally_stable2",
#     "Bergman_unstable_baseline2", "Bergman_sharp_switch_unstable5"])

#     local (; expt_id, paleorun, t_ts, O_ts, P_ts, dPdt_line, dOdt_line, P_grid, O_grid, eqb_point, fold_points, mopb, ANOX) =
#         P_O_columns_FigS1["modern_"*exptroot]
#     O2_U = PB.get_data(paleorun.output, "ocean.O2_U")
#     # O2_U_accu = cumsum(O2_U[end]) ./ collect(range(start=1, stop=100, length=100))
#     Plots.plot!(pp, accumulated_Corg_burial_norm, O2_U[end], label="PO "*expt_id, ylims=(0, 1.0), 
#                 xlabel="Accumulated_Corg_burial_norm", ylabel="O2 Utilization", color=colors[j], linestyle=linestyles[j], linewidth=linewidths[j])
# end

# #### bottom mid Corg weighted anoxia
# p_anoxia=Plots.plot()
# P_x = collect(range(0.0, 8.0, length=100))
# for (j, exptroot) in enumerate(["Bergman_constCPsea1", "Bergman_constCPsea2", "Bergman_lowO2U", "Bergman_lowCPsea", "Bergman_marginally_stable1", "Bergman_marginally_stable2",
#     "Bergman_unstable_baseline2", "Bergman_sharp_switch_unstable5"])
#     for (i, prefix) in enumerate(["modern_"]) # , "preCambrian_", "preCambrian_Corgbf0p5_"
#         local (; expt_id, paleorun, t_ts, O_ts, P_ts, dPdt_line, dOdt_line, P_grid, O_grid, eqb_point, fold_points, mopb, ANOX) =
#         P_O_columns_FigS1[prefix*exptroot]
#         Plots.plot!(p_anoxia, P_x, ANOX, color=colors[j], linestyle=linestyles[j], linewidth=linewidths[j],
#             xlabel="P_norm", ylabel="anoxia_frac", label=false)
#     end
# end

# #### bottom right the Pburial vs. P_norm
# ppp = plot()
# P_x = collect(range(0.0, 8.0, length=100))

# for (j, exptroot) in enumerate(["Bergman_constCPsea1", "Bergman_constCPsea2", "Bergman_lowO2U", "Bergman_lowCPsea", "Bergman_marginally_stable1", "Bergman_marginally_stable2",
#     "Bergman_unstable_baseline2", "Bergman_sharp_switch_unstable5"])

#     local (; expt_id, paleorun, t_ts, O_ts, P_ts, dPdt_line, dOdt_line, P_grid, O_grid, eqb_point, fold_points, mopb, ANOX) =
#         P_O_columns_FigS1["modern_"*exptroot]
#     Plots.plot!(ppp, P_x, mopb, xlims=(0,6), ylims=(0,1.6e11), color=colors[j], linestyle=linestyles[j], linewidth=linewidths[j],
#         # color=:green, linestyle=linestyles[i], label="O_n = $(1.0)", 
#         label="PO "*expt_id, xlabel="P_norm", ylabel="P burial (mol/yr)")
# end

#### plot phase planes with 3x O nullclines #######
function plot_phase_plane(exptroot)
    p=Plots.plot()

    # high CorgP x 3 nullclines
    # modern_Bergman preCambrian_Bergman preCambrian_Bergman_Corgbf0p5
    (; expt_id, t_ts, O_ts, P_ts, dPdt_line, dOdt_line, eqb_point, fold_points) = P_O_columns_FigS1["preCambrian_$(exptroot)"]
    Pbal_P, Pbal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_line)
    Obal_P, Obal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_line)
    p = Plots.plot(Pbal_P, Pbal_O, xlims=[0.0, 5.0], ylims=[0.0, 2.25], label=false, xlabel="P_norm", ylabel="O (PAL)");

    for (i, prefix) in enumerate(["modern_", "preCambrian_", "preCambrian_Corgbf0p5_"])
        local (; expt_id, paleorun, t_ts, O_ts, P_ts, dPdt_line, dOdt_line, P_grid, O_grid, eqb_point, fold_points, mopb, ANOX) =
            P_O_columns_FigS1[prefix*exptroot]
        Pbal_P, Pbal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_line)
        Obal_P, Obal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_line)
        Plots.plot!(p, Obal_P, Obal_O, color=:red, linestyle=linestyles[i], label="PO "*expt_id)

        @info "eqb_point = $(eqb_point)"
        #(; element_counts, start_point_index, end_point_index, sign_change) = SolverFunctionsOOEOAE2.find_periodic(P_ts, O_ts, t_ts; Spec_P=P_ts[end])
        if (end_point_index > 0) && (start_point_index > 0) # if the case is stable you may find the end_point_index==NaN
            Plots.plot!(p, P_ts[1:end_point_index], O_ts[1:end_point_index], color=:green, linestyle=linestyles[i], label=false)
        else
           Plots.plot!(p, P_ts, O_ts, color=:green, linestyle=linestyles[i], label=false)
        end
        if !isempty(eqb_point)
            Plots.scatter!(p, [only(eqb_point)[1]], [only(eqb_point)[2]]; color=:red, markershape=:circle, markerstrokewidth=0, label=nothing);
        end
        Plots.scatter!(p,  map(x->x[1], fold_points), map(x->x[2], fold_points); color=:blue, markershape=:circle, markerstrokewidth=0, label=nothing);
    end

    return p
end


# plot phase plane with single O nullcline
function plot_phase_plane_single_O(exptroot; ylims=[0.0, 1.5])
    p_single_O=Plots.plot()
        
    (; expt_id, t_ts, O_ts, P_ts, dPdt_line, dOdt_line, eqb_point, fold_points) = P_O_columns_FigS1[exptroot]

    @info "eqb_point = $(eqb_point)"
    (; element_counts, start_point_index, end_point_index, sign_change) = SolverFunctionsOOEOAE2.find_periodic(P_ts, O_ts, t_ts, eqb_point[1];)

    Pbal_P, Pbal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_line)
    Obal_P, Obal_O = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_line)
    p_single_O = Plots.plot(Pbal_P, Pbal_O; xlims=[0.0, 5.0], ylims, label="PO "*expt_id, xlabel="P_norm", ylabel="O (PAL)")
    Plots.plot!(p_single_O, Obal_P, Obal_O, color=:red, linestyle=:solid, label=false)
    Plots.plot!(p_single_O, P_ts, O_ts, color=:green, linestyle=:solid, label=false)
    Plots.scatter!(p_single_O, [only(eqb_point)[1]], [only(eqb_point)[2]]; color=:red, markershape=:circle, markerstrokewidth=0, label=nothing)
    Plots.scatter!(p_single_O,  map(x->x[1], fold_points), map(x->x[2], fold_points); color=:blue, markershape=:circle, markerstrokewidth=0, label=nothing);

    return p_single_O
end


#### time series, The chosen P_O_norm in single one ####
function plot_time_series(exptroot)
    p_time_series=Plots.plot()
    for (i, prefix) in enumerate(["preCambrian_"])
        local (; t_ts, O_ts, P_ts,) = P_O_columns_FigS1[prefix*exptroot]

        Plots.plot!(p_time_series, t_ts, P_ts, xlims=(0, 2.5e7), color=:blue, xlabel="tmodel (yr)", 
            linestyle=linestyles[2], linewidth=linewidths[2], label="P_norm")
        # O_norm
        Plots.plot!(p_time_series, t_ts, O_ts, xlims=(0, 2.5e7), ylims=(0,4), color=:red, xlabel="tmodel (yr)", 
            linestyle=linestyles[2], linewidth=linewidths[2], label="O_norm")
    end
    return p_time_series
end

##### time series for the waveform ######
function plot_time_series_waveform(
    exptroot;
    xlims=(0.0, 5.0e7),
    linestyle=:solid,
)
    p_time_series=Plots.plot()
    local (; t_ts, O_ts, P_ts, ) = P_O_columns_FigS1[exptroot]

    Plots.plot!(p_time_series, t_ts, P_ts; xlims, color=:blue, xlabel="tmodel (yr)", 
        linestyle, label="P_norm")
    # O_norm
    Plots.plot!(p_time_series, t_ts, O_ts; xlims, ylims=(0,4), color=:red, xlabel="tmodel (yr)", 
        linestyle, label="O_norm")
    return p_time_series
end

# https://discourse.julialang.org/t/blank-subplot-in-plots-jl/16453/2
plot_blank = Plots.plot(; legend=false, grid=false, foreground_color_subplot=:white)

p_CPsea=Plots.plot()
l = @layout[
#            [grid(4,2)];
    [grid(2, 2)];
]
p_CPsea = Plots.plot(
    plot_phase_plane("Bergman_constCPsea1"),    plot_phase_plane("Bergman_constCPsea2"), 
#     plot_time_series("Bergman_constCPsea1"),    plot_time_series("Bergman_constCPsea2"), 
    plot_phase_plane("Bergman_lowO2U"),         plot_phase_plane("Bergman_lowCPsea"),
#    plot_time_series("Bergman_lowO2U"),         plot_time_series("Bergman_lowCPsea"); 
    layout = l, size=(800, 800), left_margin = 2Plots.mm, right_margin = 2Plots.mm, bottom_margin = 1Plots.mm,
    # , left_margin = 5Plots.mm, bottom_margin = 5Plots.mm,
)
Plots.savefig(p_CPsea,  joinpath(output_figures_dir, "PO_secular_stability_SI_CPsea_20231202.svg"))

p_waveform=Plots.plot()
l = @layout[
            [grid(3,2)];
]
p_waveform = Plots.plot(
    plot_phase_plane_single_O("modern_Bergman_sharp_switch_unstable5_waveform"), 
        plot_time_series_waveform("modern_Bergman_sharp_switch_unstable5_waveform"), 
    plot_phase_plane_single_O("preCambrian_Bergman_sharp_switch_unstable5_waveform"), 
        plot_time_series_waveform("preCambrian_Bergman_sharp_switch_unstable5_waveform"),
    plot_phase_plane_single_O("preCambrian_Corgbf0p5_Bergman_sharp_switch_unstable5_waveform"),
        plot_time_series_waveform("preCambrian_Corgbf0p5_Bergman_sharp_switch_unstable5_waveform"),
    layout = l, size=(800, 600), left_margin = 2Plots.mm, right_margin = 2Plots.mm, bottom_margin = 1Plots.mm,
    # , left_margin = 5Plots.mm, bottom_margin = 5Plots.mm
)
Plots.savefig(p_waveform,  joinpath(output_figures_dir, "PO_secular_stability_SI_waveform_20231202.svg"))

p_fix_intersection=Plots.plot()
l = @layout[
            [grid(6,2)];
]
xlims_fi = (0, 2.5e7)
p_fix_intersection = Plots.plot(
    plot_phase_plane_single_O("preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness0"),
        plot_time_series_waveform("preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness0"; xlims=xlims_fi), 
    plot_phase_plane_single_O("preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness4"),
        plot_time_series_waveform("preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness4"; xlims=xlims_fi),
    plot_phase_plane_single_O("preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness3"),
        plot_time_series_waveform("preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness3"; xlims=xlims_fi), 
    plot_phase_plane_single_O("preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness2"),
        plot_time_series_waveform("preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness2"; xlims=xlims_fi),
    plot_phase_plane_single_O("preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness1"),
        plot_time_series_waveform("preCambrian_Bergman_sharp_switch_unstable5_waveform_sharpness1"; xlims=xlims_fi),
    plot_phase_plane_single_O("preCambrian_Bergman_sharp_switch_unstable5_waveform"),    
        plot_time_series_waveform("preCambrian_Bergman_sharp_switch_unstable5_waveform"; xlims=xlims_fi), 
    layout = l, 
    size=(800, 1200), 
    left_margin = 4Plots.mm, right_margin = 2Plots.mm, bottom_margin = 1Plots.mm,
    # , left_margin = 5Plots.mm, bottom_margin = 5Plots.mm,
)
Plots.savefig(p_fix_intersection,  joinpath(output_figures_dir, "PO_secular_stability_SI_fix_intersection_20231202.svg"))

p_fix_intersection_lowCPsea=Plots.plot()
l = @layout[
            [grid(7,2)];
]
xlims_fi = (0, 2.5e7)
p_fix_intersection_lowCPsea = Plots.plot(
    plot_phase_plane_single_O("preCambrian_Bergman_lowCPsea"),
        plot_time_series_waveform("preCambrian_Bergman_lowCPsea"; xlims=xlims_fi), 
    plot_phase_plane_single_O("preCambrian_Bergman_lowCPsea_sharpness1"),
        plot_time_series_waveform("preCambrian_Bergman_lowCPsea_sharpness1"; xlims=xlims_fi),
    plot_phase_plane_single_O("preCambrian_Bergman_lowCPsea_sharpness2"),
        plot_time_series_waveform("preCambrian_Bergman_lowCPsea_sharpness2"; xlims=xlims_fi), 
    plot_phase_plane_single_O("preCambrian_Bergman_lowCPsea_sharpness3"),
        plot_time_series_waveform("preCambrian_Bergman_lowCPsea_sharpness3"; xlims=xlims_fi),
    plot_phase_plane_single_O("preCambrian_Bergman_lowCPsea_sharpness4"),
        plot_time_series_waveform("preCambrian_Bergman_lowCPsea_sharpness4"; xlims=xlims_fi),
        plot_phase_plane_single_O("preCambrian_Bergman_lowCPsea_sharpness5"),
        plot_time_series_waveform("preCambrian_Bergman_lowCPsea_sharpness5"; xlims=xlims_fi),
    plot_phase_plane_single_O("preCambrian_Bergman_lowCPsea_sharpness6"),
        plot_time_series_waveform("preCambrian_Bergman_lowCPsea_sharpness6"; xlims=xlims_fi),
    layout = l, 
    size=(800, 1200), 
    left_margin = 4Plots.mm, right_margin = 2Plots.mm, bottom_margin = 1Plots.mm,
    # , left_margin = 5Plots.mm, bottom_margin = 5Plots.mm,
)
Plots.savefig(p_fix_intersection_lowCPsea,  joinpath(output_figures_dir, "PO_secular_stability_SI_fix_intersection_lowCPsea_20231217.svg"))

