using Logging
import DataFrames
import Interpolations

using Plots, Interpolations
import GLMakie # import not using so doesn't conflict with Plots
# import MarchingCubes
# import GeometryBasics
# using Roots

import PALEOboxes as PB
import PALEOmodel
import PALEOcopse

# using Interpolations

# ReactionsOOEOAE_dev add Cisotopes, link DIC_sms -= local_Corgburial!
# include("../../src/EPOC_reactions.jl")
include("ReactionsOOEOAE_dev.jl") 
include("../../src/CarbBurial_dev.jl") 
include("../../src/Uranium.jl")

include("../../src/SolverFunctionsOOEOAE2.jl")

include("ooeoae_expts.jl")

# Local figures
isdir("figures") || mkdir("figures")
dropbox_output_dir = "figures"

output_figures_dir = joinpath(dropbox_output_dir, "Major_run")
isdir(output_figures_dir) || mkdir(output_figures_dir)

######################################
# key parameters that control plot appearance
include("expt_plot3D_endP.jl")
P_lims = (0.0, 3.0)
A_lims=(1.0, 2.2)
O_lims=(0.0, 3.0)
##########################################

#####################################################
# simplified dA/dt = oxidw + ocdeg + carbw + ccdeg - mocb - mccb
# -locb = 0, no -sfw
# modified after ../OOEOAE_ZHL/6_OOEOAE_columns_copsereload_Pw_test
#####################################################

model = PB.create_model_from_config(
    joinpath(@__DIR__, "P_O_A_U_columns.yaml"), 
    "model_endP", 
    modelpars=Dict("CGconstant"=>false), # , "Aconstant"=>true
)

###############################################
# Set experiment (parameters)
###############################################
# basic setup 
ooeoae_expts(
    model, [
        # ("land_flux_Bergman", false),
        ("k_O2_U", (0.125, 1.0)), # globally stable
        ("corg_burial_fac", 1.0), # mid of two fold lines
        ("CO2pulse", [0, 5e7, 7.0e7, 1e8], [0, 0, 0, 0]), 
        ("k_anox", 100.0),
    ]
)

# Experiments as listed in Table2 in README
expts_table = [

    ("POAU 1402", "earlyT_Psilw_Uplift_VEG_stable", 
        [
            ("set_forcing_EndP"),
            # ("set_restoring_EndP"),
            ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
            ("COPSE_locb", 1380.0, 4.5e12, 0.1), 
            ("CPU_U", 40e6, 0.42*6e6, 0.42*34e6),
            # ("COPSE_landw", 2.3275e12, 4.3225e12, 13.35e12, 7.75e12, 39e9, 0.15), # original value
            #               basw   granw  carbw    oxidw   k10_phosw k_15
            ("COPSE_landw", 3e12, 9e12, 13.35e12, 7.75e12, 39e9, 0.15), 
            ("COPSE_sedcrust", 1.0, -24.0, 5.09e12, 1.25e12),
            #                  ccdeg  ocdeg

            # ("POAinit", (1.3091542834245113, 1.0710637720864389, 1.25)), 
            ("CPsea", 350.0, 350.0*4),
            ("k_O2_U", (0.99, 1.0)),
            ("corg_burial_fac", 1.3), 
            ("Uplift_gradually_recover", -252e6, -244e6, 1.1, 1.1),
            ("set_interp_forcing", "land", "force_VEG", [-1e30, -252e6 - 1, -252e6, -248e6, -246e6, -237e6], [1.0, 1.0, 0.0, 0.0, 0.2, 0.2],),
            ("set_interp_forcing", "ocean", "force_kO2U", [-1e30, -252e6 - 1, -252e6, -248e6, -248e6+1, -237e6], [0.1, 0.1, 0.1, 0.1, 0.1, 0.1],),
        ],
    ),


    ("POAU 1405", "earlyT_Psilw_Uplift_VEG_unstable", 
        [
            ("set_forcing_EndP"),
            # ("set_restoring_EndP"),
            ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
            ("COPSE_locb", 1380.0, 4.5e12, 0.1), 
            ("CPU_U", 40e6, 0.42*6e6, 0.42*34e6),
            # ("COPSE_landw", 2.3275e12, 4.3225e12, 13.35e12, 7.75e12, 39e9, 0.15), # original value
            #               basw   granw  carbw    oxidw   k10_phosw k_15
            ("COPSE_landw", 3e12, 9e12, 13.35e12, 7.75e12, 39e9, 0.15), 
            ("COPSE_sedcrust", 1.0, -24.0, 5.09e12, 1.25e12),
            #                  ccdeg  ocdeg

            # ("POAinit", (1.3091542834245113, 1.0710637720864389, 1.25)), 
            ("CPsea", 350.0, 350.0*4),
            ("k_O2_U", (0.99, 1.0)),
            ("corg_burial_fac", 1.3), 
            ("Uplift_gradually_recover", -252e6, -244e6, 1.1, 1.1),
            ("set_interp_forcing", "land", "force_VEG", [-1e30, -252e6 - 1, -252e6, -248e6, -246e6, -237e6], [1.0, 1.0, 0.0, 0.0, 0.2, 0.2],),
            ("set_interp_forcing", "ocean", "force_kO2U", [-1e30, -252e6 - 1, -252e6, -248e6, -248e6+1, -237e6], [0.1, 0.1, 0.45, 0.45, 0.1, 0.1],),
        ],
    ),
    
]

tstart = -254e6
tend =  -244e6
P_O_A_U_columns_endP = Dict()

function plot_3D_Fig4(poand, P_lims, O_lims, A_lims, A_ts, O_ts, P_ts, t_ts, poand_jac_P_P, fileroot, stage)
    # get nullcline surfaces
    (;dPdt_surf, dOdt_surf, dAdt_surf, P_grid, O_grid, A_grid) = SolverFunctionsOOEOAE2.find_nullclines_POA(
        poand; P_lims, O_lims, A_lims,
    )

    # dPdt surface
    dPdt_surf_P, dPdt_surf_O, dPdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_surf)
 
    # dA/dt = 0 surface
    dAdt_surf_P, dAdt_surf_O, dAdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dAdt_surf)
    # dA/dt = 0 intersection with dP/dt = 0   
    # dAdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dPdt_surf))
    dAdt_dPdt_line_POA = SolverFunctionsOOEOAE2.join_lines_nan_sep(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dPdt_surf))
    # split at folds
    dAdt_dPdt_linesegs_POA = SolverFunctionsOOEOAE2.split_line_sign_f(dAdt_dPdt_line_POA, poand_jac_P_P)
 
    # dO/dt = 0 surface
    dOdt_surf_P, dOdt_surf_O, dOdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_surf)
    # dO/dt = 0 intersection with dP/dt = 0
    # dOdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[2], dPdt_surf))
    dOdt_dPdt_line_POA = SolverFunctionsOOEOAE2.join_lines_nan_sep(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[2], dPdt_surf))
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

    ########################################
    # make 3D phase plane
    # find where dA/dt=0 along line with dO/dt=0 and dP/dt=0
    eqb_point = only(SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dOdt_dPdt_line_POA; verbose=true))
    @info "eqb_point POA: $eqb_point"

    # P, O phase plane at constant A: find dP/dt = 0 for a representative value of A
    # (can get this from surfaces we have already calculated)
    A_target_val = eqb_point[3] # roughly midpoint of limit cycle
    A_ks = findfirst(x -> x > A_target_val, A_grid)
    A_val = A_grid[A_ks] # actual A used
    dPdt_linesegs_const_A = SolverFunctionsOOEOAE2.split_line_sign_f(dPdt_surf[:, A_ks], poand_jac_P_P)
    Obal_P, Obal_O, _ = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_surf[:, 1]) # dO/dt = 0 is independent of A so just pick first A

    # P, A plane at constant O: find dP/dt = 0, dA/dt = 0 for one value of O
    # (need to calculate nullclines in P, A surface at constant O)
    # O_val = 1.5 # roughly midpoint of limit cycle
    O_val = eqb_point[2] # O value of the eqb point!
    dPdt_line_const_O, dAdt_line_const_O = SolverFunctionsOOEOAE2.find_nullclines_const_O(poand, O_val; P_lims, A_lims)
    # split at folds
    dPdt_linesegs_const_O = SolverFunctionsOOEOAE2.split_line_sign_f(dPdt_line_const_O, poand_jac_P_P)
    
    # O, A plane, "projecting" critical manifold ie using P value(s) from critical manifold at each O, A
    dAdt_dPdt_linesegs_POA_end = Vector{Any}()
   
    # make 3D phase plane
    fig, axs, pltobj = plot_3D_critical_manifold(dPdt_surf, A_grid, folds_dPdt_lines_POA)

    # # dA/dt = 0 surface
    # plot_3D_A_nullcline(axs, dAdt_surf)
    SolverFunctionsOOEOAE2.plot_segments!(
        GLMakie.lines!, axs, dAdt_dPdt_linesegs_POA, [1, 2, 3], (:dash, :solid), nothing;
        color=:black, linewidth=glm_linewidth,
    )

    # # dO/dt = 0 surface
    # plot_3D_O_nullcline(axs, dOdt_surf)
    SolverFunctionsOOEOAE2.plot_segments!(
        GLMakie.lines!, axs, dOdt_dPdt_linesegs_POA, [1, 2, 3], (:dash, :solid), nothing;
        color=:red, linewidth=glm_linewidth,
    )

    # TIME series
    # GLMakie.lines!(axs, P_ts, O_ts, A_ts; color=:green, linewidth=1*glm_linewidth,)
    if stage == 1
        # plot from 254 to 252 in red
        indices = findall(x -> -254e6 ≤ x ≤ -252e6, t_ts)
        GLMakie.lines!(axs, P_ts[indices], O_ts[indices], A_ts[indices]; color=:red, linewidth=2.5*glm_linewidth,)
    elseif stage == 2
        # plot from 254 to 252 in red, 252-248 in blue
        indices = findall(x -> -254e6 ≤ x ≤ -252e6, t_ts)
        GLMakie.lines!(axs, P_ts[indices], O_ts[indices], A_ts[indices]; color=:red, linewidth=2.5*glm_linewidth,)
        indices = findall(x -> -252e6 ≤ x ≤ -248e6, t_ts)
        GLMakie.lines!(axs, P_ts[indices], O_ts[indices], A_ts[indices]; color=:green, linewidth=2.5*glm_linewidth,)
    elseif stage == 3
        # plot from 254 to 252 in red, 252-248 in blue, 248-244 in green
        indices = findall(x -> -254e6 ≤ x ≤ -252e6, t_ts)
        GLMakie.lines!(axs, P_ts[indices], O_ts[indices], A_ts[indices]; color=:red, linewidth=2.5*glm_linewidth,)
        indices = findall(x -> -252e6 ≤ x ≤ -248e6, t_ts)
        GLMakie.lines!(axs, P_ts[indices], O_ts[indices], A_ts[indices]; color=:green, linewidth=2.5*glm_linewidth,)
        indices = findall(x -> -248e6 ≤ x ≤ -244e6, t_ts)
        GLMakie.lines!(axs, P_ts[indices], O_ts[indices], A_ts[indices]; color=:blue, linewidth=2.5*glm_linewidth,)
    end
    GLMakie.scatter!(axs, eqb_point; color=:red, markersize=glm_eqb_markersize,)

    # display(fig)
    GLMakie.save(joinpath(output_figures_dir, "POAU_"*fileroot*".png"), fig; px_per_unit=5,)
end


# we need the phase plane for 253, 250, 245 Ma
for (expt_id, fileroot, vector_pars) in expts_table

# (expt_id, fileroot, vector_pars) = expts_table[1]

    ooeoae_expts(
        model, vector_pars
    )

    tspan = (tstart, tend) # yr # tspan=(-1000e6, 0)

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
    global poand_253 = SolverFunctionsOOEOAE2.POAnormDeriv(model, modeldata; include_jacobian, fixed_values=initial_state, fixed_t=-253e6,)
    global poand_250 = SolverFunctionsOOEOAE2.POAnormDeriv(model, modeldata; include_jacobian, fixed_values=initial_state, fixed_t=-250e6,)
    global poand_245 = SolverFunctionsOOEOAE2.POAnormDeriv(model, modeldata; include_jacobian, fixed_values=initial_state, fixed_t=-245e6,)
    if include_jacobian
        # d/dP_m(dP_n/dt) component of Jacobian to identify folds
        # -ve is an attracting surface, +ve is repelling
        global poand_jac_P_P_253(P_n, O_n, A_n) = SolverFunctionsOOEOAE2.jacobian(poand_253, P_n, O_n, A_n)[1, 1]
        global poand_jac_P_P_250(P_n, O_n, A_n) = SolverFunctionsOOEOAE2.jacobian(poand_250, P_n, O_n, A_n)[1, 1]
        global poand_jac_P_P_245(P_n, O_n, A_n) = SolverFunctionsOOEOAE2.jacobian(poand_245, P_n, O_n, A_n)[1, 1]
    else
        # label everything as 'attracting'
        global poand_jac_P_P_253(P_n, O_n, A_n) = -1.0 
        global poand_jac_P_P_250(P_n, O_n, A_n) = -1.0 
        global poand_jac_P_P_245(P_n, O_n, A_n) = -1.0 
    end

    # get nullcline surfaces
    plot_3D_Fig4(poand_253, P_lims, O_lims, A_lims, A_ts, O_ts, P_ts, t_ts, poand_jac_P_P_253, fileroot*"_253", 1)
    plot_3D_Fig4(poand_250, P_lims, O_lims, A_lims, A_ts, O_ts, P_ts, t_ts, poand_jac_P_P_250, fileroot*"_250", 2)
    plot_3D_Fig4(poand_245, P_lims, O_lims, A_lims, A_ts, O_ts, P_ts, t_ts, poand_jac_P_P_245, fileroot*"_245", 3)

    ##########################################
    # save everything in 3D for the sum plot!
    P_O_A_U_columns_endP[fileroot] = 
        (; expt_id, paleorun, modeldata, t_ts, A_ts, P_ts, O_ts)
    ##########################################
end # end of for loop


# # # ####### For the 2D TS part #######
function plot_time_series_POA(exptname)
    local (; paleorun, t_ts, A_ts, P_ts, O_ts) = 
        P_O_A_U_columns_endP[exptname]

    # @info "eqb_point POA: $eqb_point"
    # @info "end point = $(P_ts[end]), $(O_ts[end]), $(A_ts[end])"

    # for time series
    pp_norm = Plots.plot()
    Plots.plot!(pp_norm, t_ts, O_ts, label="O_norm", linestyle=linestyles[1], color=:red)
    Plots.plot!(pp_norm, t_ts, P_ts, label="P_norm", linestyle=linestyles[1], color=:blue)
    Plots.plot!(pp_norm, t_ts, A_ts, label="A_norm", linestyle=linestyles[1], color=:black)

    d13C = PB.get_data(paleorun.output, "ocean.DIC_delta")

    p_d13C = Plots.plot(t_ts, d13C, label="ocean.DIC_delta", linestyle=linestyles[1], color=:black, ylims=(-2,8)) # , xlims=(-254e6, -244e6)

    l = @layout[
        grid(2,1);
    ]

    p_sum = Plots.plot(pp_norm, p_d13C,
                    layout = l, left_margin = 5Plots.mm, right_margin = 5Plots.mm, bottom_margin = 5Plots.mm, size=(1000, 300)) # 

    savefig(p_sum, joinpath(output_figures_dir, "POAU_"*exptname*"time_series_20240905.svg"))
    return pp_norm, p_d13C # 
end

plot_time_series_POA("earlyT_Psilw_Uplift_VEG_stable")
plot_time_series_POA("earlyT_Psilw_Uplift_VEG_unstable")
