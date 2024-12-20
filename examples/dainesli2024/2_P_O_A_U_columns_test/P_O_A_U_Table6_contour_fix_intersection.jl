using Logging
import DataFrames
# import Interpolations
import LinearAlgebra

using Plots
import GLMakie # import not using so doesn't conflict with Plots
# import MarchingCubesP
# import GeometryBasics
# using Roots

import PALEOboxes as PB
import PALEOmodel
# import PALEOreactions
import PALEOcopse

# ReactionsOOEOAE_dev add Cisotopes, link DIC_sms -= local_Corgburial!
include("../../../src/EPOC_reactions.jl")

include("../../../src/SolverFunctionsOOEOAE2.jl")

include("../ooeoae_expts.jl")


# Archived figures
# dropbox_output_dir = joinpath(@__DIR__, "../../figures/2_P_O_A_U_columns_test")

# Local figures
isdir("figures") || mkdir("figures")
dropbox_output_dir = "figures"

output_figures_dir = joinpath(dropbox_output_dir, "P_O_A_U_Table6_contour_20231212")
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
    modelpars=Dict("CGconstant"=>true), # , "Aconstant"=>true
)

###############################################
# Set experiment (parameters)
###############################################
# basic setup

CP_sea_ox, CP_sea_anox = 115.4, 461.6
ooeoae_expts(
    model, [
        # ("land_flux_Bergman", false),
        ("k_O2_U", (0.23, 0.25)), # sharp switch
        # ("CPsea", 115.4, 461.6), # base setup values for Bergman_COPSE
        ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
        ("CO2pulse", [0, 5e7, 7.0e7, 1e8], [0, 0, 0, 0]), 
        ("COPSE_locb", 1000.0, 0.0),
        ("k_anox", 100.0),
        ("OPAinit", (0.548*3.7e19, 1.21*3.1e15, 4.53*3.193e18)),
        ("CPsea", CP_sea_ox, CP_sea_anox),
    ]
)

# See  P-O case PALEOexamples/src/OOEOAE_base/1_P_O_columns_test/P_O_columns_FigS1_20231202.jl
# 
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
calc_k_O2_U1(k_O2_U_2) = (0.561 - 0.529*k_O2_U_2)/(4*0.471)

# fold location (corg_burial_fac) calculated using
# PALEOexamples/src/OOEOAE_base/2_P_O_A_U_columns_test/P_O_A_U_test_folds_nullcline.jl

fold_sharpness_grid = [
    # k_O2_U        # fold corg_burial_frac 
    ((0.23, 0.25),   (1.1379837770625496, 0.4779556059928637)), # sharp switch case
    ((0.2065, 0.325), (1.226048, 0.460539)),
    ((0.185, 0.4),   (1.3188879002408516, 0.4422216626440729)),
    ((0.157, 0.5),   (1.42257014899951, 0.41866169623822364)),
    # ((0.145, 0.53),  (1.4332971063051407, 0.4123448129006402)),  # typo ?
    ((0.149, 0.53),  (1.441415, 0.414498)),
    ((0.561/4, 0.561), (1.49, 0.4)),  # no fold ! (this is the neutrally stable case)
    ((0.132, 0.59),  (1.53, 0.39)),   # no fold !
]

# calculate "sharpness factor" (to use as y axis of plot) and get fold locations
y_Sharpness, extfold_1_cbf, extfold_2_cbf = [], [], [] # fold + extension into stable regime
fold_sharpness, fold_1_cbf, fold_2_cbf = [], [], []  # only where sharpness >= 0
for fold_pars in fold_sharpness_grid
    (k_O2_U, f_cbf) = fold_pars
    sharpness = SolverFunctionsOOEOAE2.find_sharpness(k_O2_U, CP_sea_anox/CP_sea_ox)
    push!(y_Sharpness, sharpness)
    push!(extfold_1_cbf, f_cbf[2]) # oxic fold + extension into stable regime
    push!(extfold_2_cbf, f_cbf[1]) # anoxic fold + extension into stable regime
    if sharpness >= 0.0 # don't include globally stable cases
        push!(fold_sharpness, sharpness)
        push!(fold_1_cbf, f_cbf[2]) # oxic fold
        push!(fold_2_cbf, f_cbf[1]) # anoxic fold
    end
end

# set corg grid relative to fold location
# corg_burial_fac_rel_grid = [0.05, 0.1, 0.33, 0.67, 0.9, 0.95] 
corg_burial_fac_rel_grid = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95] 
# corg_burial_fac_rel_grid = [0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975] 

# allocate array to store output (oscillation period)
# NB: the x grid is relative to the fold location
z_Periodicity_rel = zeros(length(fold_sharpness_grid), length(corg_burial_fac_rel_grid)) # numerical period 
jac_rel = Matrix{Any}(undef, size(z_Periodicity_rel)...) # to hold Jacobian matrices
linear_growth_rate_rel = fill(NaN, size(jac_rel)) #  growth rate from jacobian (-ve for damped oscillations)
linear_period_rel = fill(NaN, size(jac_rel)) # period from linear stability analysis


for (i, fold_pars) in enumerate(fold_sharpness_grid)
    (k_O2_U, f_cbf) = fold_pars

    for (j, corg_burial_fac_rel) in enumerate(corg_burial_fac_rel_grid)
        # convert corg relative to fold to an actual value
        # NB f_cbf[2] is the oxic fold (low corg_burial_fac)
        corg_burial_fac = f_cbf[2] + corg_burial_fac_rel*(f_cbf[1] - f_cbf[2])

        ooeoae_expts(
            model, [("k_O2_U", k_O2_U), ("corg_burial_fac", corg_burial_fac),]
        )

        tspan = (0.0, 1e8) # yr # tspan=(-1000e6, 0)

        #########################################################
        # Initialize
        #########################################################

        initial_state, modeldata = PALEOmodel.initialize!(model)

        ################################
        # get time_series
        ####################################
        (paleorun, A_ts, O_ts, P_ts) = SolverFunctionsOOEOAE2.find_time_series(
            model, initial_state, modeldata;
            has_A=true, tspan
        )
        t_ts = PB.get_data(paleorun.output, "global.tforce")


        poand = SolverFunctionsOOEOAE2.POAnormDeriv(model, modeldata; include_jacobian=true)

        eqb_point = SolverFunctionsOOEOAE2.find_intersection_point_newton(poand)
        @info "eqb_point Newton: $eqb_point"

        # Jacobian for linear stability analysis
        jac_rel[i, j] = SolverFunctionsOOEOAE2.jacobian(poand, eqb_point...)
        ev = LinearAlgebra.eigvals(jac_rel[i, j])
        for e in ev
            # we expect a pair of complex-valued eigenvalues,
            # real part gives growth rate, 2π / (imag part) gives oscillation period
            if imag(e) != 0
                linear_growth_rate_rel[i, j] = real(e)
                linear_period_rel[i, j] = 2*π/abs(imag(e)) # only valid for damped oscillations
                break
            end
        end
        
        @info "eqb_point jacobian: " jac_rel[i, j]
        @info "eqb_point linear growth rate: $(linear_growth_rate_rel[i, j]) (yr) period $(linear_period_rel) (yr)"
       
        # find the period numerically
        (; element_counts, start_point_index, end_point_index, sign_change) = 
            SolverFunctionsOOEOAE2.find_periodic(
                P_ts, O_ts, t_ts, eqb_point;
                # Spec_P=eqb_point[1],
                omit_t_start=20e6, # omit initial transient when estimating P to use for zero-crossing
                # plot=true, # for debugging
            )

        @info "find_periodic returned element_counts ", element_counts
        if !isempty(element_counts)
            periodic = argmax(element_counts)
        else
            @warn "element_counts is empty, setting periodic to NaN"
            periodic = NaN
        end
        z_Periodicity_rel[i, j] = periodic
    end

end # end of for loop

# Summary of expts to overlay on regime
# TODO move this to a common location
expt_summary = [
    # id        k_O2_U          corg_burial_fac
    ("6",       (0.05, 0.75),   0.8), # stable
    ("7",       (0.23, 0.25),   1.1),
    ("8",       (0.23, 0.25),   0.7),
    ("9",       (0.23, 0.25),   0.5),
    ("10",      (0.23, 0.25),   1.25), # OAE, also S35, S36, S37
    ("12",      (0.23, 0.25),   0.45), # OOE
    ("S33",     (0.125, 0.25),   0.7), # fast limit cycle NB: this is not actually comparable ! (it isnt in the "fix intersection point" grid)
]
expt_id, expt_sharpness, expt_corg = [], [], []
for (id, k_O2_U, corg) in expt_summary
    push!(expt_id, id)
    push!(expt_sharpness, SolverFunctionsOOEOAE2.find_sharpness(k_O2_U, CP_sea_anox/CP_sea_ox))
    push!(expt_corg, corg)
end

gr(size=(500, 500)) # we require an 8.7cm wide single-column figure 
levels = [2.0, 2.5, 3.0, 5.0, 7.5, 10.0]

# period from numerical model compared to linear stability analysis
p_fig6_rel = contour(
    corg_burial_fac_rel_grid, reverse(y_Sharpness), reverse(z_Periodicity_rel./1e6; dims=1);
    cbar=false, clabels=true, color=:black, levels=levels, label="numerical",
);
# period from linear stability analysis
contour!(
    p_fig6_rel,
    corg_burial_fac_rel_grid, reverse(y_Sharpness), reverse(linear_period_rel./1e6; dims=1);
    cbar=false, clabels=true, color=:black, linestyle=:dash, levels=levels, label="linear",
);
# plot zero-level contour of growth rate as a crude way to use linear interpolation to find where it should be zero
contour!(
    p_fig6_rel, corg_burial_fac_rel_grid, reverse(y_Sharpness), reverse(linear_growth_rate_rel; dims=1);
    cbar=false, clabels=false, color=:blue, linestyle=:dash, levels=[0.0],
)
plot!(
    p_fig6_rel;  
    xlims=(0.0, 1.0), xlabel="corg_burial_fac relative to folds", 
    ylabel="Sharpness",
)
display(p_fig6_rel)
savefig(p_fig6_rel, joinpath(output_figures_dir, "Fig6_Stability_diagram_Corg_rel.svg"))

# 2D interpolation to map relative corg_burial into absolute values
corg_burial_fac_grid, y_new_Sharpness, z_Periodicity = 
    SolverFunctionsOOEOAE2.resample_map_grid(
        corg_burial_fac_rel_grid, y_Sharpness, z_Periodicity_rel, extfold_1_cbf, extfold_2_cbf;
    )

# growth rate from jacobian
# 2D interpolation to map relative corg_burial into absolute values
_, _, linear_growth_rate = 
    SolverFunctionsOOEOAE2.resample_map_grid(
        corg_burial_fac_rel_grid, y_Sharpness, linear_growth_rate_rel, extfold_1_cbf, extfold_2_cbf;
    )

p_fig6 = contour(
    corg_burial_fac_grid, reverse(y_new_Sharpness), reverse(z_Periodicity./1e6; dims=1);
    xlims=(0.0, 1.8), xlabel="carbon burial efficiency (O nullcline)", xticks=[0.5, 1.0, 1.5],
    ylims=(-0.3, 1.0), ylabel="fold sharpness (P nullcline)", yticks=[-0.25, 0.0, 0.25, 0.5, 0.75, 1.0],
    cbar=false, clabels=true, color=:black, levels=levels,
);
# plot zero-level contour of growth rate as a crude way to use linear interpolation to find where it should be zero
contour!(
    p_fig6, corg_burial_fac_grid, reverse(y_new_Sharpness), reverse(linear_growth_rate; dims=1);
    cbar=false, clabels=false, color=:blue, linestyle=:dash, levels=[0.0],
)

plot!(p_fig6, fold_1_cbf, fold_sharpness; color=:brown, label=nothing);
plot!(p_fig6, fold_2_cbf, fold_sharpness; color=:green, label=nothing);
plot!(p_fig6, [last(fold_1_cbf), last(fold_2_cbf)], last(fold_sharpness).*[1, 1]; color=:blue, label=nothing);
scatter!(p_fig6, expt_corg, expt_sharpness; color=:red, label=nothing)

display(p_fig6)
savefig(p_fig6, joinpath(output_figures_dir, "Fig6_Stability_diagram_Corg.svg"))