using Logging
import DataFrames

using Plots
import XLSX
import Interpolations

using Plots, Interpolations
import PALEOboxes as PB
import PALEOmodel
import PALEOcopse
import PALEOocean

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

# ######################################
# # key parameters that control plot appearance
# include("expt_plot3D_endP.jl")
# ##########################################

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
    ("POAU 1401", "earlyT_Psilw_Uplift_VEG_stable_min", 
        [
            ("set_forcing_EndP_0915"),
            # ("set_restoring_EndP"),
            ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
            ("COPSE_locb", 1000.0, 4.5e12, 0.1),
            ("CPsea", 180.0, 180.0*4), 
            ("CPU_U", 40e6, 0.76*6e6, 0.76*34e6),
            # ("COPSE_landw", 2.3275e12, 4.3225e12, 13.35e12, 7.75e12, 39e9, 0.15), # original value
            #               basw   granw  carbw    oxidw   k10_phosw k_15
            ("COPSE_landw", 1.5e12, 4.5e12, 13.35e12, 7.75e12, 39e9, 0.15), 
            ("COPSE_sedcrust", 1.0, -22.0, 5.0e12, 1.25e12),
            #                  ccdeg  ocdeg

            # ("POAinit", (1.3091542834245113, 1.0710637720864389, 1.25)), 
            ("k_O2_U", (0.99, 1.0)),
            ("corg_burial_fac", 1.0), 
            ("Uplift_gradually_recover", -252e6, -244e6, 0.8, 0.8),
            ("set_interp_forcing", "land", "force_VEG", [-1e30, -252e6 - 1, -252e6, -248e6, -246e6, -237e6], [1.0, 1.0, 0.0, 0.0, 0.25, 0.25],),
            ("set_interp_forcing", "ocean", "force_kO2U", [-1e30, -252e6 - 1, -252e6, -248e6, -248e6+1, -237e6], [0.1, 0.1, 0.01, 0.01, 0.1, 0.1],),
        ],
    ),

    ("POAU 1402", "earlyT_Psilw_Uplift_VEG_stable", 
        [
            ("set_forcing_EndP_0915"),
            # ("set_restoring_EndP"),
            ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
            ("COPSE_locb", 1000.0, 4.5e12, 0.1),
            ("CPsea", 180.0, 180.0*4), 
            ("CPU_U", 40e6, 0.76*6e6, 0.76*34e6),
            # ("COPSE_landw", 2.3275e12, 4.3225e12, 13.35e12, 7.75e12, 39e9, 0.15), # original value
            #               basw   granw  carbw    oxidw   k10_phosw k_15
            ("COPSE_landw", 1.5e12, 4.5e12, 13.35e12, 7.75e12, 39e9, 0.15), 
            ("COPSE_sedcrust", 1.0, -22.0, 5.0e12, 1.25e12),
            #                  ccdeg  ocdeg

            # ("POAinit", (1.3091542834245113, 1.0710637720864389, 1.25)), 
            ("k_O2_U", (0.99, 1.0)),
            ("corg_burial_fac", 1.0), 
            ("Uplift_gradually_recover", -252e6, -244e6, 0.8, 0.8),
            ("set_interp_forcing", "land", "force_VEG", [-1e30, -252e6 - 1, -252e6, -248e6, -246e6, -237e6], [1.0, 1.0, 0.0, 0.0, 0.25, 0.25],),
            ("set_interp_forcing", "ocean", "force_kO2U", [-1e30, -252e6 - 1, -252e6, -248e6, -248e6+1, -237e6], [0.1, 0.1, 0.1, 0.1, 0.1, 0.1],),
        ],
    ),

    ("POAU 1403", "earlyT_Psilw_Uplift_VEG_stable_max", 
        [
            ("set_forcing_EndP_0915"),
            # ("set_restoring_EndP"),
            ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
            ("COPSE_locb", 1000.0, 4.5e12, 0.1),
            ("CPsea", 180.0, 180.0*4), 
            ("CPU_U", 40e6, 0.76*6e6, 0.76*34e6),
            # ("COPSE_landw", 2.3275e12, 4.3225e12, 13.35e12, 7.75e12, 39e9, 0.15), # original value
            #               basw   granw  carbw    oxidw   k10_phosw k_15
            ("COPSE_landw", 1.5e12, 4.5e12, 13.35e12, 7.75e12, 39e9, 0.15), 
            ("COPSE_sedcrust", 1.0, -22.0, 5.0e12, 1.25e12),
            #                  ccdeg  ocdeg

            # ("POAinit", (1.3091542834245113, 1.0710637720864389, 1.25)), 
            ("k_O2_U", (0.99, 1.0)),
            ("corg_burial_fac", 1.0), 
            ("Uplift_gradually_recover", -252e6, -244e6, 0.8, 0.8),
            ("set_interp_forcing", "land", "force_VEG", [-1e30, -252e6 - 1, -252e6, -248e6, -246e6, -237e6], [1.0, 1.0, 0.0, 0.0, 0.25, 0.25],),
            ("set_interp_forcing", "ocean", "force_kO2U", [-1e30, -252e6 - 1, -252e6, -248e6, -248e6+1, -237e6], [0.1, 0.1, 0.2, 0.2, 0.1, 0.1],),
        ],
    ),


    ("POAU 1404", "earlyT_Psilw_Uplift_VEG_unstable_min", 
        [
            ("set_forcing_EndP_0915"),
            # ("set_restoring_EndP"),
            ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
            ("COPSE_locb", 1000.0, 4.5e12, 0.1),
            ("CPsea", 180.0, 180.0*4), 
            ("CPU_U", 40e6, 0.76*6e6, 0.76*34e6),
            # ("COPSE_landw", 2.3275e12, 4.3225e12, 13.35e12, 7.75e12, 39e9, 0.15), # original value
            #               basw   granw  carbw    oxidw   k10_phosw k_15
            ("COPSE_landw", 1.5e12, 4.5e12, 13.35e12, 7.75e12, 39e9, 0.15), 
            ("COPSE_sedcrust", 1.0, -22.0, 5.0e12, 1.25e12),
            #                  ccdeg  ocdeg

            # ("POAinit", (1.3091542834245113, 1.0710637720864389, 1.25)), 
            ("k_O2_U", (0.99, 1.0)),
            ("corg_burial_fac", 1.0), 
            ("Uplift_gradually_recover", -252e6, -244e6, 0.8, 0.8),
            ("set_interp_forcing", "land", "force_VEG", [-1e30, -252e6 - 1, -252e6, -248e6, -246e6, -237e6], [1.0, 1.0, 0.0, 0.0, 0.25, 0.25],),
            ("set_interp_forcing", "ocean", "force_kO2U", [-1e30, -252e6 - 1, -252e6, -248e6, -248e6+1, -237e6], [0.1, 0.1, 0.3, 0.3, 0.1, 0.1],),
        ],
    ),

    ("POAU 1405", "earlyT_Psilw_Uplift_VEG_unstable", 
        [
            ("set_forcing_EndP_0915"),
            # ("set_restoring_EndP"),
            ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
            ("COPSE_locb", 1000.0, 4.5e12, 0.1),
            ("CPsea", 180.0, 180.0*4), 
            ("CPU_U", 40e6, 0.76*6e6, 0.76*34e6),
            # ("COPSE_landw", 2.3275e12, 4.3225e12, 13.35e12, 7.75e12, 39e9, 0.15), # original value
            #               basw   granw  carbw    oxidw   k10_phosw k_15
            ("COPSE_landw", 1.5e12, 4.5e12, 13.35e12, 7.75e12, 39e9, 0.15), 
            ("COPSE_sedcrust", 1.0, -22.0, 5.0e12, 1.25e12),
            #                  ccdeg  ocdeg

            # ("POAinit", (1.3091542834245113, 1.0710637720864389, 1.25)), 
            ("k_O2_U", (0.99, 1.0)),
            ("corg_burial_fac", 1.0), 
            ("Uplift_gradually_recover", -252e6, -244e6, 0.8, 0.8),
            ("set_interp_forcing", "land", "force_VEG", [-1e30, -252e6 - 1, -252e6, -248e6, -246e6, -237e6], [1.0, 1.0, 0.0, 0.0, 0.25, 0.25],),
            ("set_interp_forcing", "ocean", "force_kO2U", [-1e30, -252e6 - 1, -252e6, -248e6, -248e6+1, -237e6], [0.1, 0.1, 0.4, 0.4, 0.1, 0.1],),
        ],
    ),

    ("POAU 1406", "earlyT_Psilw_Uplift_VEG_unstable_max", 
        [
            ("set_forcing_EndP_0915"),
            # ("set_restoring_EndP"),
            ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
            ("COPSE_locb", 1000.0, 4.5e12, 0.1),
            ("CPsea", 180.0, 180.0*4), 
            ("CPU_U", 40e6, 0.76*6e6, 0.76*34e6),
            # ("COPSE_landw", 2.3275e12, 4.3225e12, 13.35e12, 7.75e12, 39e9, 0.15), # original value
            #               basw   granw  carbw    oxidw   k10_phosw k_15
            ("COPSE_landw", 1.5e12, 4.5e12, 13.35e12, 7.75e12, 39e9, 0.15), 
            ("COPSE_sedcrust", 1.0, -22.0, 5.0e12, 1.25e12),
            #                  ccdeg  ocdeg

            # ("POAinit", (1.3091542834245113, 1.0710637720864389, 1.25)), 
            ("k_O2_U", (0.99, 1.0)),
            ("corg_burial_fac", 1.0), 
            ("Uplift_gradually_recover", -252e6, -244e6, 0.8, 0.8),
            ("set_interp_forcing", "land", "force_VEG", [-1e30, -252e6 - 1, -252e6, -248e6, -246e6, -237e6], [1.0, 1.0, 0.0, 0.0, 0.25, 0.25],),
            ("set_interp_forcing", "ocean", "force_kO2U", [-1e30, -252e6 - 1, -252e6, -248e6, -248e6+1, -237e6], [0.1, 0.1, 0.5, 0.5, 0.1, 0.1],),
        ],
    ),
    
]

isdir(output_figures_dir) || mkdir(output_figures_dir)
tstart = -254e6
tend =  -244e6

P_O_A_U_columns_earlyT = Dict()

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

    ##########################################
    # save everything in 3D for the sum plot!
    P_O_A_U_columns_earlyT[fileroot] = 
        (; expt_id, paleorun, t_ts)
    ##########################################

end # end of for loop

linestyles=[:dashdotdot, :solid, :dash, :dashdotdot, :solid, :dash]
linewidths=[0.1, 2, 0.1, 0.1, 2, 0.1]
colors=[:blue, :blue, :blue, :orangered3, :orangered3, :orangered3]

plot_forcing=plot(); plot_stability=plot(); plot_anox=plot(); plot_TEMP=(); plot_C=plot(); plot_U=plot()
(plot_PTB_C, plot_PTB_U, plot_PTB_TEMP, df_C, df_U, df_TEMP) = plot_data_C_U_O("PTB-C-U-O.xlsx")

plot_PTB_C = scatter(df_C.Age/-1e6, df_C.d13Ccarb, marker=:circle,color=:grey100, alpha=0.01, ylabel="δ¹³Ccarb", label=false,)
plot_PTB_U = scatter(df_U.Age/-1e6, df_U.d238U, marker=:circle, color=:gray, alpha=0.05, ylabel = "δ²³⁸Ucarb", label=false,)

for (i, fileroot) in enumerate(["earlyT_Psilw_Uplift_VEG_stable_min", "earlyT_Psilw_Uplift_VEG_stable", "earlyT_Psilw_Uplift_VEG_stable_max",
                                    "earlyT_Psilw_Uplift_VEG_unstable_min", "earlyT_Psilw_Uplift_VEG_unstable", "earlyT_Psilw_Uplift_VEG_unstable_max"])
    (; expt_id, paleorun, t_ts) = P_O_A_U_columns_earlyT[fileroot]
    VEG = PB.get_data(paleorun.output, "land.VEG")
    UPLIFT = PB.get_data(paleorun.output, "global.UPLIFT")
    O2_U = PB.get_data(paleorun.output, "ocean.kO2U")
    TEMP = PB.get_data(paleorun.output, "global.TEMP")
    d13C = PB.get_data(paleorun.output, "ocean.DIC_delta")
    d238U = PB.get_data(paleorun.output, "ocean.U_delta")
    anox = PB.get_data(paleorun.output, "oceanfloor.anoxia_burial_frac")

    t_ts = t_ts /-1e6

    plot!(plot_forcing, t_ts, VEG, color=:green, xlabel="", label=false, linewidth=2)
    plot!(plot_forcing, t_ts, UPLIFT, color=:orange, grid = false, xlabel="", ylabel="Forcings", xlims=(244, 254), xflip=true, label=false, linewidth=2)

    plot!(plot_stability, t_ts, O2_U, grid = false, color=colors[i], linewidth=linewidths[i], linestyle=linestyles[i], ylabel="Stability", label=false, xlims=(244, 254), ylims=(0,0.5), xflip=true,)

    plot!(plot_PTB_C, t_ts, d13C, grid = false, xlabel="", ylabel="δ¹³Ccarb", color=colors[i], linewidth=linewidths[i], linestyle=linestyles[i], label=false, xlims=(244, 254), xflip=true, ylims=(-6, 10))

    plot!(plot_PTB_U, t_ts, d238U, grid = false, xlabel="tmodel (Ma)", ylabel="δ²³⁸Ucarb", color=colors[i], linewidth=linewidths[i], linestyle=linestyles[i], label=false, ylims=(-1.2, 0.4), xlims=(244, 254), xflip=true,)
    
    plot!(plot_anox, t_ts, anox, grid = false, xlabel="tmodel (Ma)", ylabel="Ocean anoxia", color=colors[i], linewidth=linewidths[i], linestyle=linestyles[i], label=false, xlims=(244, 254), xflip=true,)
end

# fill the range
function plot_range!(pp, fileroot, data; color, for_temp=false)
    (; paleorun, t_ts) = P_O_A_U_columns_earlyT[fileroot*"_min"]
    line_min = PB.get_data(paleorun.output, data)
    t_ts_min = t_ts

    (; paleorun, t_ts) = P_O_A_U_columns_earlyT[fileroot*"_max"]
    line_max = PB.get_data(paleorun.output, data)
    t_ts_max = t_ts

    itp_min = interpolate((t_ts_min,), line_min, Gridded(Linear()))
    itp_max = interpolate((t_ts_max,), line_max, Gridded(Linear()))

    ts = -254e6:1e5:-244e6

    min_interp = itp_min.(ts)
    max_interp = itp_max.(ts)

    ts = ts /-1e6

    if for_temp
        plot!(twinx(), ts, min_interp, fillrange=max_interp, fillalpha=0.1, fillcolor=color, linecolor=:transparent, label=false, xlims=(244, 254), xflip=true, ylims=(285.0, 297.0), yflip=false)
    else
        plot!(pp, ts, min_interp, fillrange=max_interp, fillalpha=0.1, fillcolor=color, linecolor=:transparent, label=false,)
    end
end

plot_range!(plot_stability, "earlyT_Psilw_Uplift_VEG_stable", "ocean.kO2U"; color=:blue)
plot_range!(plot_stability, "earlyT_Psilw_Uplift_VEG_unstable", "ocean.kO2U"; color=:orangered3)

plot_range!(plot_PTB_C, "earlyT_Psilw_Uplift_VEG_stable", "ocean.DIC_delta"; color=:blue)
plot_range!(plot_PTB_C, "earlyT_Psilw_Uplift_VEG_unstable", "ocean.DIC_delta"; color=:orangered3)

plot_range!(plot_PTB_U, "earlyT_Psilw_Uplift_VEG_stable", "ocean.U_delta"; color=:blue)
plot_range!(plot_PTB_U, "earlyT_Psilw_Uplift_VEG_unstable", "ocean.U_delta"; color=:orangered3)

plot_range!(plot_anox, "earlyT_Psilw_Uplift_VEG_stable", "oceanfloor.anoxia_burial_frac"; color=:blue)
plot_range!(plot_anox, "earlyT_Psilw_Uplift_VEG_unstable", "oceanfloor.anoxia_burial_frac"; color=:orangered3)



# the TEMP and d18O plot has double-yaxis, need to handle seperately 
plot_TEMP = scatter(df_TEMP.Age/-1e6, df_TEMP.d18O, grid = false, marker=:circle, color=:gray, yguidefontcolor=:gray, alpha=0.05, ylabel = "δ¹⁸Oapatite", label=false, ylims=(17, 21), yflip=true, xlims=(244, 254), xflip=true,)
for (i, fileroot) in enumerate(["earlyT_Psilw_Uplift_VEG_stable_min", "earlyT_Psilw_Uplift_VEG_stable", "earlyT_Psilw_Uplift_VEG_stable_max",
                                    "earlyT_Psilw_Uplift_VEG_unstable_min", "earlyT_Psilw_Uplift_VEG_unstable", "earlyT_Psilw_Uplift_VEG_unstable_max"])
    (; expt_id, paleorun, t_ts) = P_O_A_U_columns_earlyT[fileroot]
    TEMP = PB.get_data(paleorun.output, "global.TEMP")
    t_ts = t_ts /-1e6
    plot!(twinx(), t_ts, TEMP, grid = false, ylabel="TEMP (k)", color=colors[i], linewidth=linewidths[i], linestyle=linestyles[i], yguidefontcolor=:orangered3, label=false, xlims=(244, 254), xflip=true, ylims=(285.0, 297.0), yflip=false)
end

plot_range!(plot_TEMP, "earlyT_Psilw_Uplift_VEG_stable", "global.TEMP"; color=:blue, for_temp=true)
plot_range!(plot_TEMP, "earlyT_Psilw_Uplift_VEG_unstable", "global.TEMP"; color=:orangered3, for_temp=true)



# sum plots
l = @layout[
            grid(3,2);
    ]

p_sum_timeseries = Plots.plot(plot_forcing, plot_stability, plot_PTB_C, plot_TEMP, plot_PTB_U, plot_anox, 
                        layout = l, left_margin = 5Plots.mm, right_margin = 5Plots.mm, size=(800, 800))

display(p_sum_timeseries)
# savefig(p_sum_timeseries, joinpath(output_figures_dir, "POAU_TS_only_20241224.svg"))  