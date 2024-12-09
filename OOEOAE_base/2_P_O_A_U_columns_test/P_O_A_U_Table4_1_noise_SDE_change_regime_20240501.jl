using Logging
import DataFrames
import Interpolations
# import DifferentialEquations
import StochasticDiffEq
import SciMLBase
import CSV
# import GLMakie

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOcopse

# using Interpolations

# ReactionsOOEOAE_dev add Cisotopes, link DIC_sms -= local_Corgburial!
include("../ReactionsOOEOAE_dev.jl")
include("../SolverFunctionsOOEOAE2.jl")

include("CarbBurial_dev.jl")

include("../ooeoae_expts.jl")
include("../ooeoae_plots.jl")

######################################
# key parameters that control plot appearance
# include("expt_plot3D.jl")

# Archived figures
dropbox_output_dir = joinpath(@__DIR__, "../../figures/2_P_O_A_U_columns_test")

# Local figures
isdir("figures") || mkdir("figures")
dropbox_output_dir = "figures"

output_figures_dir = joinpath(dropbox_output_dir, "P_O_A_U_Table7_white_noise_change_regimes")
isdir(output_figures_dir) || mkdir(output_figures_dir)

#####################################################
# simplified dA/dt = oxidw + ocdeg + carbw + ccdeg - mocb - mccb
# -locb = 0, no -sfw
# modified after ../OOEOAE_ZHL/6_OOEOAE_columns_copsereload_Pw_test
#####################################################

model = PB.create_model_from_config(
    joinpath(@__DIR__, "P_O_A_U_columns.yaml"), 
    "model3", 
    modelpars=Dict(
        "CGconstant"=>false, 
        #  "Aconstant"=>true,
        "state_norm"=>true,
        # "state_norm"=>false,
    )
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
    ]
)

# Experiments as listed in Table2 in README
expts_table = [

    # across the oxic threshold
    ("POA 2010", "preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_CO2_whitenoise", [("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
        # ("CO2pulse", age, pulse), 
        ("OPAinit", (0.4*3.7e19, 0.9*3.1e15, 4.07*3.193e18)),
        ("corg_burial_fac", 1.0), # Corgb oxic fold ~= 1.0, Corgb anoxic fold ~= 0.5
        ("set_interp_forcing", "ocean", "force_corg_bf",    [0.0, 1e8, 1e8+1, 2e8, 2e8+1, 3e8, 3e8+1, 4e8, 4e8+1, 5e8], [0.9, 0.9, 1.1, 1.1, 1.12, 1.12, 1.2, 1.2, 1.4, 1.4]),
    ]),

    # # across the anoxic threshold
    # ("POA 2011", "preCambrian_Bergman_sharp_switch_unstable5_Corgb5_Psilw_only_CO2pulse_whitenoise", [("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
    #     # ("CO2pulse", age, pulse), 
    #     ("OPAinit", (0.8369220346325966*3.7e19, 3.7066508042082518*3.1e15, 4.072859528806692*3.193e18)),
    #     ("corg_burial_fac", 1.0), # Corgb oxic fold ~= 1.0, Corgb anoxic fold ~= 0.5
    #     ("set_interp_forcing", "ocean", "force_corg_bf",    [0.0, 1e8, 1e8+1, 2e8, 2e8+1, 3e8, 3e8+1, 4e8, 4e8+1, 5e8], [0.4, 0.4, 0.45, 0.45, 0.485, 0.485, 0.5, 0.5, 0.7, 0.7]),
    # ]),

]

isdir(output_figures_dir) || mkdir(output_figures_dir)

P_O_A_U_columns_Table4_1 = Dict()

for (expt_id, fileroot, vector_pars) in expts_table

# (expt_id, fileroot, vector_pars) = expts_table[1]

    ooeoae_expts(
        model, vector_pars
    )

    tspan = (0.0, 5e8) # yr # tspan=(-1000e6, 0)

    #########################################################
    # Initialize
    #########################################################

    initial_state, modeldata = PALEOmodel.initialize!(model)

    # call ODE function to check derivative
    initial_deriv = similar(initial_state)
    PALEOmodel.ODE.ModelODE(modeldata)(initial_deriv, initial_state , nothing, 0.0)
    println("initial_state", initial_state)
    println("initial_deriv", initial_deriv)
    println("integrate, ODE")
    ################################
    # get time_series
    ####################################
    # (paleorun, A_ts, O_ts, P_ts) = SolverFunctionsOOEOAE2.find_time_series(model, initial_state, modeldata; has_A=true, tspan, dtmax=1e4)
    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    m = PALEOmodel.ODE.ModelODE(modeldata; solver_view=modeldata.solver_view_all, dispatchlists=modeldata.dispatchlists_all)
    # f = SciMLBase.ODEFunction{true}(m)

    # # using ODEProblem
    # reltol = 1e-5
    # alg = Sundials.CVODE_BDF()
    # prob = SciMLBase.ODEProblem(m, initial_state, tspan, nothing)
    # @time sol = SciMLBase.solve(prob, alg; reltol=reltol)
    # PALEOmodel.ODE.print_sol_stats(sol)
    # paleorun.output = PALEOmodel.OutputWriters.OutputMemory()
    # PALEOmodel.ODE.calc_output_sol!(paleorun.output, paleorun.model, sol, tspan, initial_state, modeldata)

    # get first index (mol total) of atmocean A (A_solve is the normalized value provided to solver)
    # (an isotope variable is represented by two values, first index is total mol, second index is mol*delta)
    A_index::Int = first(PB.get_indices(modeldata.solver_view_all.stateexplicit, "atmocean.A_solve"))

    amplitude_noise1 = [0.0, 0.0, 0.0]/sqrt(1e6)
    # amplitude_noise2 = [0.05, 0.05, 0.05]/sqrt(1e6)
    amplitude_noise2 = [0.1, 0.1, 0.1]/sqrt(1e6)
    amplitude_noise3 = [0.2, 0.2, 0.2]/sqrt(1e6)
    amplitude_noise4 = [0.3, 0.3, 0.3]/sqrt(1e6)

    amplitude_noise_v = [amplitude_noise1, amplitude_noise2, amplitude_noise3, amplitude_noise4]

    for i in [1]
        amplitude_noise = amplitude_noise_v[i]

        function g(du, u, p, t)
            # NB: we are now using normalized A_solve, so scale noise appropriately
            # assume characteristic time for A reservoir is ~1e6 yr, then "standard noise" random perturbation will be ~sqrt(1e6) after 1e6 yr
            # so du[A_index] = 1/sqrt(1e6) will be ~comparable to a perturbation of A by A_norm (or A_solve by 1.0)
            if t <= tspan[end]/3
                du[A_index] = amplitude_noise[1] # small noise
            elseif tspan[end]/3 < t <= tspan[end] * 2/3
                du[A_index] = amplitude_noise[2] # mid noise
            elseif tspan[end] * 2/3 < t <= tspan[end]
                du[A_index] = amplitude_noise[3] # large noise
            end
            # @info "$du"
        end
        # a = SciMLBase.SDEFunction{true}(m, g)

        # alg = DifferentialEquations.SOSRA()
        alg = StochasticDiffEq.SOSRA()
        prob = SciMLBase.SDEProblem(m, g, initial_state, tspan, nothing)
        @time sol = SciMLBase.solve(prob, alg, save_noise=true) # , reltol=1e-9， dt = 1e4
        PALEOmodel.ODE.print_sol_stats(sol)
        paleorun.output = PALEOmodel.OutputWriters.OutputMemory()
        PALEOmodel.ODE.calc_output_sol!(paleorun.output, paleorun.model, sol, tspan, initial_state, modeldata)
        
        # sol.W is the Wiener process, which is intergrated by white noise
        noise_A = []
        for i in 1:length(sol)
            push!(noise_A, sol.W[i][A_index])
        end
        noise_A_diff = diff(noise_A) # the white noise

        for i in 1:length(sol.t)-1
            t = sol.t[i]
            if t <= tspan[end]/3
                noise_A_diff[i] = noise_A_diff[i] * amplitude_noise[1] # small noise
            elseif tspan[end]/3 < t <= tspan[end] * 2/3
                noise_A_diff[i] = noise_A_diff[i] * amplitude_noise[2] # mid noise
            elseif tspan[end] * 2/3 < t <= tspan[end]
                noise_A_diff[i] = noise_A_diff[i] * amplitude_noise[3] # large noise
            end
        end

        l = @layout[
            grid(6,1)
        ]

        p_sum = Plots.plot(
            plot(sol.t[1:end-1], noise_A_diff, label=false, ylabel="Noise", xlabel="",),
            plot(paleorun.output, "ocean.corg_bf", label=false, xlabel="", ylims=(0.25, 1.4)),
            plot(paleorun.output, "atmocean.A_norm", label=false, xlabel="",),
            plot(paleorun.output, "atmocean.O_norm", label=false, xlabel="",),
            plot(paleorun.output, "ocean.P_norm", label=false, xlabel="",),
            plot(paleorun.output, "ocean.DIC_delta", label=false),
            layout = l, 
            left_margin = 5Plots.mm, bottom_margin = 1Plots.mm,
            size=(400, 800)
        )
        t_ts = PB.get_data(paleorun.output, "global.tforce")
        P_ts = PB.get_data(paleorun.output, "ocean.P_norm")
        O_ts = PB.get_data(paleorun.output, "atmocean.O_norm")
        A_ts = PB.get_data(paleorun.output, "atmocean.A_norm")
        DIC_delta = PB.get_data(paleorun.output, "ocean.DIC_delta")

        display(p_sum)
        savefig(p_sum, joinpath(output_figures_dir, fileroot * "_$(amplitude_noise[1])" * ".svg"))

        
        # ###########################################################
        # # phase plane plot with nullclines and folds
        # ########################################################
        # include_jacobian = true
        # # poand: initial_state, poand_end: end_state
        # poand = SolverFunctionsOOEOAE2.POAstatenormDeriv(model, modeldata; fixed_values=initial_state, fixed_t=tspan[1], include_jacobian)
        # if include_jacobian
        #     # d/dP_m(dP_n/dt) component of Jacobian to identify folds
        #     # -ve is an attracting surface, +ve is repelling
        #     global poand_jac_P_P(P_n, O_n, A_n) = SolverFunctionsOOEOAE2.jacobian(poand, P_n, O_n, A_n)[1, 1]
        #     # global poand_jac_P_P_end(P_n, O_n, A_n) = SolverFunctionsOOEOAE2.jacobian(poand_end, P_n, O_n, A_n)[1, 1]
        # else
        #     # label everything as 'attracting'
        #     global poand_jac_P_P(P_n, O_n, A_n) = -1.0 
        #     # global poand_jac_P_P_end(P_n, O_n, A_n) = -1.0 
        # end

        # # get nullcline surfaces
        # (;dPdt_surf, dOdt_surf, dAdt_surf, A_grid) = SolverFunctionsOOEOAE2.find_nullclines_POA(poand; P_lims, O_lims, A_lims=(2.0, A_lims[2]))
    
    
        # ########################################################
        # # dA/dt = 0 surface, for initial_state
        # dAdt_surf_P, dAdt_surf_O, dAdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dAdt_surf)
        # # dA/dt = 0 intersection with dP/dt = 0   
        # dAdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dPdt_surf))
        # # split at folds
        # dAdt_dPdt_linesegs_POA = SolverFunctionsOOEOAE2.split_line_sign_f(dAdt_dPdt_line_POA, poand_jac_P_P)

        # ########################################################

        # ########################################################
        # # dO/dt = 0 surface
        # dOdt_surf_P, dOdt_surf_O, dOdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_surf)
        # # dO/dt = 0 intersection with dP/dt = 0
        # dOdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[2], dPdt_surf))
        # # split at folds
        # dOdt_dPdt_linesegs_POA = SolverFunctionsOOEOAE2.split_line_sign_f(dOdt_dPdt_line_POA, poand_jac_P_P)
        # ########################################################
        
        # # find folds (P component of Jacobian = 0)
        # folds_dPdt_lines_POA = SolverFunctionsOOEOAE2.Isoline.find_isolines(poand_jac_P_P, dPdt_surf)
        # # order folds so first index is lowest P
        # if length(folds_dPdt_lines_POA) == 2
        #     f_1, f_2 = folds_dPdt_lines_POA[1], folds_dPdt_lines_POA[2]
        #     P_1, P_2 = first(f_1)[1], first(f_2)[1]
        #     if P_1 > P_2
        #         folds_dPdt_lines_POA[1:2] .= f_2, f_1
        #     end
        # end

        # # find where dA/dt=0 along line with dO/dt=0 and dP/dt=0
        # eqb_point = only(SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dOdt_dPdt_line_POA; verbose=true))
        # @info "eqb_point POA: $eqb_point"

        # A_target_val = eqb_point[3] # roughly midpoint of limit cycle
        # A_ks = findfirst(x -> x > A_target_val, A_grid)
        # A_val = A_grid[A_ks] # actual A used
        # dPdt_linesegs_const_A = SolverFunctionsOOEOAE2.split_line_sign_f(dPdt_surf[:, A_ks], poand_jac_P_P)
        # Obal_P, Obal_O, _ = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_surf[:, 1]) # dO/dt = 0 is independent of A so just pick first A

        # O_val = eqb_point[2] # roughly midpoint of limit cycle
        # dPdt_line_const_O, dAdt_line_const_O = SolverFunctionsOOEOAE2.find_nullclines_const_O(poand, O_val; P_lims, A_lims)
        # # split at folds
        # dPdt_linesegs_const_O = SolverFunctionsOOEOAE2.split_line_sign_f(dPdt_line_const_O, poand_jac_P_P)

        # fig, axs, pltobj = plot_3D_critical_manifold(dPdt_surf, A_grid, folds_dPdt_lines_POA)

        # # add intersection lines of O_nullcline and slow-manifold surface
        # SolverFunctionsOOEOAE2.plot_segments!(
        #     GLMakie.lines!, axs, dOdt_dPdt_linesegs_POA, [1, 2, 3], (:dash, :solid), nothing; 
        #     color=:red, linewidth=glm_linewidth,
        # )
        # # add intersection lines of A_nullcline and slow-manifold surface
        # SolverFunctionsOOEOAE2.plot_segments!(
        #     GLMakie.lines!, axs, dAdt_dPdt_linesegs_POA, [1, 2, 3], (:dash, :solid), nothing;
        #     color=:black, linewidth=glm_linewidth,
        # )
        # # add limit cycle time series
        # (; element_counts, start_point_index, end_point_index) = SolverFunctionsOOEOAE2.find_periodic(P_ts, O_ts, t_ts, eqb_point) 
        # GLMakie.lines!(
        #     axs,
        #     P_ts, O_ts, A_ts;
        #     # P_ts[start_point_index:end_point_index], O_ts[start_point_index:end_point_index], A_ts[start_point_index:end_point_index];
        #     color=:green, linestyle=:solid, linewidth=2*glm_linewidth,
        # )

        # GLMakie.scatter!(axs, eqb_point; color=:red, markersize=glm_eqb_markersize,)
        # GLMakie.scatter!(axs, (P_ts[end], O_ts[end], A_ts[end]); color=:blue, markersize=glm_eqb_markersize,)

        # GLMakie.save(joinpath(output_figures_dir, 
        #     fileroot * "_$(amplitude_noise[1])" * ".png"), 
        #     fig; px_per_unit=5, # Makie 0.20 increase resolution of saved figure (600 x 5 = 3000 pixels)
        # ) # closes window

        ##########################################
        # save everything in 3D for the sum plot!
        P_O_A_U_columns_Table4_1[fileroot * "_$(amplitude_noise[1])"] = 
            (; expt_id, paleorun, t_ts, A_ts, P_ts, O_ts, DIC_delta,
                # A_val, O_val,
                # dPdt_surf, A_grid, dAdt_surf, dOdt_surf,
                # dPdt_linesegs_const_A, Obal_P, Obal_O, 
                # dPdt_line_const_O, dAdt_line_const_O, dPdt_linesegs_const_O,
                # folds_dPdt_lines_POA, dOdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA, dAdt_dPdt_linesegs_POA_end, eqb_point
            )
        ##########################################
    end
end


# (;expt_id, paleorun, t_ts, A_ts, P_ts, O_ts, DIC_delta,) = P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_CO2_whitenoise_0.0"]
# plot_oxic_0 = plot(t_ts, DIC_delta, ylabel="δ¹³C (‰)", color=:orange, label=false)
# (;expt_id, paleorun, t_ts, A_ts, P_ts, O_ts, DIC_delta,) = P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb2_Psilw_only_CO2_whitenoise_0.0003"]
# corg_bf = PB.get_data(paleorun.output, "ocean.corg_bf")
# plot_oxic_corg_bf = plot(t_ts, corg_bf, label=false, xlabel="", ylims=(0.8, 1.4))
# plot_oxic_0003 = plot(t_ts, DIC_delta, ylabel="δ¹³C (‰)", color=:blue, label=false)

# (;expt_id, paleorun, t_ts, A_ts, P_ts, O_ts, DIC_delta,) = P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb5_Psilw_only_CO2pulse_whitenoise_0.0"]
# plot_anoxic_0 = plot(t_ts, DIC_delta, ylabel="δ¹³C (‰)", color=:orange, label=false)
# (;expt_id, paleorun, t_ts, A_ts, P_ts, O_ts, DIC_delta,) = P_O_A_U_columns_Table4_1["preCambrian_Bergman_sharp_switch_unstable5_Corgb5_Psilw_only_CO2pulse_whitenoise_0.0003"]
# corg_bf = PB.get_data(paleorun.output, "ocean.corg_bf")
# plot_anoxic_corg_bf = plot(t_ts, corg_bf, label=false, xlabel="", ylims=(0.2, 0.8))
# plot_anoxic_0003 = plot(t_ts, DIC_delta, ylabel="δ¹³C (‰)", color=:blue, label=false)

#         l = @layout[
#             grid(6,1)
#         ]

#         p_sum = Plots.plot(
#             plot_oxic_corg_bf,
#             plot_oxic_0, plot_oxic_0003,
#             plot_anoxic_corg_bf,
#             plot_anoxic_0, plot_anoxic_0003,
#             layout = l, 
#             left_margin = 5Plots.mm, bottom_margin = 1Plots.mm,
#             size=(400, 800)
#         )
#         display(p_sum)
#         savefig(p_sum, joinpath(output_figures_dir, "Whitenoise_oxic_anoxic_sum_20240511.svg"))