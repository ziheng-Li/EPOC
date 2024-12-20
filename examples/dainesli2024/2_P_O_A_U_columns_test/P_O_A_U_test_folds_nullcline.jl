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

output_figures_dir = joinpath(dropbox_output_dir, "P_O_A_U_Table6")
isdir(output_figures_dir) || mkdir(output_figures_dir)

######################################
# key parameters that control plot appearance
include("expt_plot3D.jl")
O_lims=(0.0, 2.0) # higher upper limit, needs to include "corner"  of dP/dt=0 surface otherwise plots fail
##########################################

#####################################################
# simplified dA/dt = oxidw + ocdeg + carbw + ccdeg - mocb - mccb
# -locb = 0, no -sfw
# modified after ../OOEOAE_ZHL/6_OOEOAE_columns_copsereload_Pw_test
#####################################################

model = PB.create_model_from_config(
    joinpath(@__DIR__, "P_O_A_U_columns.yaml"), 
    "model1", 
    # modelpars=Dict("CGconstant"=>false), # , "Aconstant"=>true
    modelpars=Dict("CGconstant"=>true),
)

###############################################
# Set experiment (parameters)
###############################################
# basic setup 
ooeoae_expts(
    model, [
        # ("k_O2_U", (0.15, 0.85)), # oscillations, not full amplitude
        # ("corg_burial_fac", 1.0), # unstable
        # ("corg_burial_fac", 1.674), # eqb point at "oxic" fold
        # ("corg_burial_fac", 0.60), # eqb point at "anoxic" fold
        ("CPsea", 115.4, 230.8), # base setup values for Bergman_COPSE
        ("OPAinit", (0.35*3.7e19, 1.13*3.1e15, 4.36*3.193e18)),
        # ("OPAinit", (0.548*3.7e19, 1.21*3.1e15, 4.53*3.193e18)),
        ("k_anox", 100.0),  # slighly less sharp transition: default 1000 causes problems with fold detection
        # ("CO2pulse", [0, 1e6, 1.01e6, 1e12], [0, 0, -7.9e12/2, -7.9e12/2]),  # halve degassing just after start
    ]
)

# One test case from PALEOexamples\src\OOEOAE_base\2_P_O_A_U_columns_test\P_O_A_U_Table2_3Dplot2.jl
# (limit cycle oscillation)
(fileroot, vector_pars) =
    # ("preCambrian_Bergman_unstable_basline1_Psilw_only",  [
    #     ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
    #     ("COPSE_locb", 1000.0, 0.0),    
    #     # ("k_O2_U", (0.30, 0.70)), # marginally unstable
    #     # ("k_O2_U", (0.23, 0.25)), # sharp switch
    #     # ("k_O2_U", (0.2065, 0.325)),
    #     # ("k_O2_U", (0.185, 0.4)), 
    #     # ("k_O2_U", (0.157, 0.5)), 
    #     # ("k_O2_U", (0.149, 0.53)),
    #     # ("k_O2_U", (0.561/4, 0.561)), 
    #     ("k_O2_U", (0.132, 0.59)),
    # ]) 
    # ("preCambrian_Bergman_sharp_switch_unstable5_waveform_Psilw_only_sharpness4", [("k_O2_U", (0.132, 0.59)), 
    #         ("corg_burial_fac", 0.65), ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), ("corg_burial_fac", 1.53), ("OPAinit", (0.548*3.7e19, 1.21*3.1e15, 4.53*3.193e18))])
    ("preCambrian_Bergman_lowCPsea",  [
        ("P_weathering", 3.9e10, 1.0, 0.0, 0.0), 
        ("COPSE_locb", 1000.0, 0.0),  
        ("corg_burial_fac", 0.7), 
        ("k_O2_U", (0.23, 0.25)),
        # ("k_O2_U", (0.193, 0.30)), 
        # ("k_O2_U", (0.179, 0.32)), 
        # ("k_O2_U", (0.336/2, 0.336)),
        # ("k_O2_U", (0.158, 0.35)), 
        # ("k_O2_U", (0.123, 0.4)),
    ]) 

ooeoae_expts(
    model, vector_pars
)

tspan = (0.0, 1e8) # yr # tspan=(-1000e6, 0)

#########################################################
# Initialize
#########################################################

default_initial_state, modeldata = PALEOmodel.initialize!(model)


###########################################################
# function object for phase plane plot with nullclines and folds
########################################################

poand = SolverFunctionsOOEOAE2.POAnormDeriv(model, modeldata; fixed_values=default_initial_state, include_jacobian=true)
# d/dP_m(dP_n/dt) component of Jacobian to identify folds
# -ve is an attracting surface, +ve is repelling
poand_jac_P_P(P_n, O_n, A_n) = SolverFunctionsOOEOAE2.jacobian(poand, P_n, O_n, A_n)[1, 1]

#########################
# modify initial state to supplied P_n, O_n, A_n

# unmodified 
initial_state = copy(default_initial_state)

# initial conditions for eqb point on anoxic fold
# initial_state = SolverFunctionsOOEOAE2.modify_state(poand, 3.95, 1.18, 4.08) 

################################
# get time_series
####################################
(paleorun, A_ts, O_ts, P_ts) = SolverFunctionsOOEOAE2.find_time_series(
    model, initial_state, modeldata; 
    has_A=true, tspan
)

# get nullcline surfaces
(; dPdt_surf, dOdt_surf, dAdt_surf, A_grid) = SolverFunctionsOOEOAE2.find_nullclines_POA(
    poand; 
    P_lims, O_lims, A_lims=(2.0, A_lims[2])
    # P_lims=(0,8), O_lims=(0,4), A_lims=(2,5)
)

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

# plot dP/dt=0
fig, axs, pltobj = plot_3D_critical_manifold(dPdt_surf, A_grid, folds_dPdt_lines_POA)


# dA/dt = 0 surface
plot_3D_A_nullcline(axs, dAdt_surf)
# dA/dt = 0 intersection with dP/dt = 0   
dAdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines(
    (P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dPdt_surf)
)
# split at folds
dAdt_dPdt_linesegs_POA = SolverFunctionsOOEOAE2.split_line_sign_f(dAdt_dPdt_line_POA, poand_jac_P_P)
SolverFunctionsOOEOAE2.plot_segments!(
    GLMakie.lines!, axs, dAdt_dPdt_linesegs_POA, [1, 2, 3], (:dash, :solid), nothing; 
    color=:black, linewidth=glm_linewidth,
)

# dO/dt = 0 surface
plot_3D_O_nullcline(axs, dOdt_surf)
# dO/dt = 0 intersection with dP/dt = 0
dOdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[2], dPdt_surf))
# split at folds
dOdt_dPdt_linesegs_POA = SolverFunctionsOOEOAE2.split_line_sign_f(dOdt_dPdt_line_POA, poand_jac_P_P)
SolverFunctionsOOEOAE2.plot_segments!(
    GLMakie.lines!, axs, dOdt_dPdt_linesegs_POA, [1, 2, 3], (:dash, :solid), nothing;
    color=:red, 
    # linewidth=2,
    linewidth=glm_linewidth,
)


# time series
GLMakie.lines!(
    axs, P_ts, O_ts, A_ts;
    # color=A_ts,
    # colormap=colormap,
    # colorrange=(2.0, 5.0),
    color=:green,
    # linewidth=4, # slow ?!
    linewidth=2*glm_linewidth,
)

# Find eqb point
# eqbpoint = only(SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dOdt_dPdt_line_POA; verbose=true))
eqbpoint = SolverFunctionsOOEOAE2.find_intersection_point_newton(poand)

###########################################################
# Find corg_burial_fac such that eqb point is on a fold
##########################################################

# Try two ways of finding where dA/dt = 0 crosses fold lines
# 1) Search along dA_n/dt=0 line in dP_n/dt =0 surface 
# (find where dA/dt=0 nullcline in dP/dt=0 surface crosses folds)
dAdt_dPdt_folds_approach_1 = SolverFunctionsOOEOAE2.Isoline.find_zeros_line(poand_jac_P_P, dAdt_dPdt_line_POA; verbose=true)
@info "dAdt_dPdt_folds_approach_1: $dAdt_dPdt_folds_approach_1"

# 2) Search along each fold line for point where dA_n/dt = 0 
dAdt_dPdt_folds_approach_2 = eltype(dAdt_dPdt_line_POA)[]
for (i, fl) in enumerate(folds_dPdt_lines_POA)
    fzl = SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], fl; verbose=true)
    if length(fzl) == 1
        push!(dAdt_dPdt_folds_approach_2, first(fzl))
    else
        @info "    dAdt_dPdt_folds_approach_2 fold $i had $(length(fzl)) crossing points for dA_n/dt=0 line"
    end
end
@info "dAdt_dPdt_folds_approach_2: $dAdt_dPdt_folds_approach_2"

GLMakie.scatter!(
    axs, SolverFunctionsOOEOAE2.tuples_to_coords(dAdt_dPdt_folds_approach_1)...;
    color=:blue,
)

display(fig)

# reset model, modeldata parameter corg_burial_fac
function reset_corg_burial_fac!(model, modeldata, corgbf_value)

    ooeoae_expts(
        model,  [("corg_burial_fac", corgbf_value)]
    )

    PB.dispatch_setup(model, :setup, modeldata)
end

# solve for corg_burial_fac in range corbf_range that gives dO/dt = 0 at dOdtzeropoint
# pond should be a function such that (dP_n/dt, dO_n/dt, dA_n/dt) = poand(P_n, O_n, A_n)
function find_corg_burial_fac!(model, modeldata, poand, corgbf_range, fold_dPdt_line)

    # find value of dO/dt at (P_n, O_n, A_n) = dOdtzeropoint, for given corg_burial_fac
    function dOdtatdAdt0(corgbf_value, p)
        reset_corg_burial_fac!(model, modeldata, corgbf_value)

        # find where dA/dt=0 nullcline crosses fold line 
        dOdtzeropoint = only(SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], fold_dPdt_line; verbose=true))
        return poand(dOdtzeropoint...)[2]
    end

    # solve for corg_burial_fac that gives dO/dt = 0 at dOdtzeropoint
    corgbf_value, niter, converged = SolverFunctionsOOEOAE2.Isoline.ITP_solve(
        dOdtatdAdt0, corgbf_range; 
        ϵ=1e-6*sum(corgbf_range)/2,
        verbose=true
    )

    converged || error("find_corg_burial_fac! failed corgbf_range $cobf_range")

    return corgbf_value
end

# find corg_burial_fac for eqb point on each fold
cobf_fold_1 = find_corg_burial_fac!(model, modeldata, poand, [0.7, 2.0], folds_dPdt_lines_POA[1])
cobf_fold_2 = find_corg_burial_fac!(model, modeldata, poand, [0.5, 0.9], folds_dPdt_lines_POA[2])

######################################################
# Find eqb points across a range of corg_burial_fac
####################################################

"""
    find_eqb_point(model, modeldata, poand, corg_burial_fac, dPdt_surf) -> (P_n, O_n, A_n)

Find eqb point for specified corg_burial_fac
    
Assumes dPdt_surf doesn't change, then recalculates O, A nullclines in that surface and finds intersection

NB: resets model to supplied corg_burial_fac parameter !!!
"""
function find_eqb_point(model, modeldata, poand, corg_burial_fac, dPdt_surf)
    reset_corg_burial_fac!(model, modeldata, corg_burial_fac)
    # Recalculate nullclines in dP/dt=0 surface
    # dO/dt = 0 intersection with dP/dt = 0
    # dOdt_dPdt_line_POA = only(SolverFunctionsOOEOAE2.Isoline.find_isolines((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[2], dPdt_surf))
    # # find where dA/dt=0 along line with dO/dt=0 and dP/dt=0
    # eqbpoint = only(SolverFunctionsOOEOAE2.Isoline.find_zeros_line((P_n, O_n, A_n)->poand(P_n, O_n, A_n)[3], dOdt_dPdt_line_POA; verbose=true))
    eqbpoint = SolverFunctionsOOEOAE2.find_intersection_point_newton(poand)

    return eqbpoint
end

# find set of eqb points and overlay
corg_burial_fac_eqb_points = collect(0.1:0.1:2.0)
stable_eqb_points = eltype(dAdt_dPdt_line_POA)[]
unstable_eqb_points = eltype(dAdt_dPdt_line_POA)[]
for corg_burial_fac in corg_burial_fac_eqb_points
    eqbpoint = find_eqb_point(model, modeldata, poand, corg_burial_fac, dPdt_surf)
    ddPdPdt = poand_jac_P_P(eqbpoint...)
    if ddPdPdt >= 0.0
        push!(unstable_eqb_points, eqbpoint)
    else
        push!(stable_eqb_points, eqbpoint)
    end
end

# https://docs.makie.org/stable/reference/plots/scatter/ for marker types
GLMakie.scatter!(
    axs, SolverFunctionsOOEOAE2.tuples_to_coords(stable_eqb_points)...;
    color=:green,
)
GLMakie.scatter!(
    axs, SolverFunctionsOOEOAE2.tuples_to_coords(unstable_eqb_points)...;
    marker=:xcross,
    color=:green,
)

display(fig)

GLMakie.save(
    joinpath(output_figures_dir, "test_OPAinit1.png"), fig;
    # joinpath(output_figures_dir, "test_OPAinit2.png"), fig;
    px_per_unit=5, # Makie 0.20 increase resolution of saved figure (600 x 5 = 3000 pixels)
)

# redisplay needed to reset scaling ?
display(fig)

@info "corg_burial_fac for eqb point on oxic fold $cobf_fold_1"
@info "corg_burial_fac for eqb point on anoxic fold $cobf_fold_2"
@info "(edit configuration to change corg_burial_fac and rerun to confirm)"