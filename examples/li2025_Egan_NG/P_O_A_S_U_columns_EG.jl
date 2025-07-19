using Logging
import DataFrames
import Interpolations

using Plots
import XLSX

import PALEOboxes as PB
import PALEOmodel
import PALEOcopse

# using Interpolations

# ReactionsOOEOAE_dev add Cisotopes, link DIC_sms -= local_Corgburial!
include("../../src/EPOC_reactions.jl")
include("ReactionsS_dev.jl")

include("../../src/SolverFunctionsOOEOAE2.jl")

include("ooeoae_expts.jl")

######################################
# # key parameters that control plot appearance
include("expt_plot_EG.jl")

# Local figures
isdir("figures") || mkdir("figures")
dropbox_output_dir = "figures"

output_figures_dir = joinpath(dropbox_output_dir, "P_O_A_U_major_run")
isdir(output_figures_dir) || mkdir(output_figures_dir)
##########################################

model = PB.create_model_from_config(
    joinpath(@__DIR__, "P_O_A_S_U_columns.yaml"), 
    # "model1",  # no S
    "model2",  # include S
    modelpars=Dict(), 
    # modelpars=Dict("restore_P"=>true, "restore_O"=>true, "restore_A"=>true, "restore_U"=>true), 
)

###############################################
# Set experiment (parameters)
###############################################
# basic setup 
ooeoae_expts(
    model, [
        # ("land_flux_Bergman", false),
        ("CPsea", 115.4, 461.6), # base setup values for Bergman_COPSE
        ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
        ("CO2pulse", [0, 5e7, 7.0e7, 1e8], [0, 0, 0, 0]), 
        ("COPSE_locb", 1000.0, 0.0),
        ("k_anox", 100.0),
    ]
)

# Experiments as listed in Table2 in README
expts_table = [
    ("POASU 1", "EG_regime_stableanoxic_to_unstable1", 
        [
            ("k_O2_U", (0.23, 0.25)), 
            ("corg_burial_fac", 1.0), 
            ("OPAinit", (0.3*3.7e19, 4*3.1e15, 4.2*3.193e18)),

            ("set_interp_forcing", "ocean", "force_corg_bf", [-1e30, 5e6 - 1, 5e6, 26e6, 26e6+1, 35e6], [0.3, 0.3, 0.7, 0.7, 1.0, 1.0],),
        ]),
    
    ("POASU 2", "EG_regime_stableanoxic_to_unstable2", 
        [
            ("k_O2_U", (0.23, 0.25)), 
            ("corg_burial_fac", 1.0), 
            ("OPAinit", (0.3*3.7e19, 4*3.1e15, 4.2*3.193e18)),

            ("set_interp_forcing", "ocean", "force_corg_bf", [-1e30, 5e6 - 1, 5e6, 26e6, 26e6+1, 35e6], [0.3, 0.3, 0.8, 0.8, 1.1, 1.1],),
        ]),

    ("POASU 3", "EG_regime_stableanoxic_to_unstable3", 
        [
            ("k_O2_U", (0.23, 0.25)), 
            ("corg_burial_fac", 1.0), 
            ("OPAinit", (0.3*3.7e19, 4*3.1e15, 4.2*3.193e18)),

            ("set_interp_forcing", "ocean", "force_corg_bf", [-1e30, 5e6 - 1, 5e6, 26e6, 26e6+1, 35e6], [0.3, 0.3, 0.6, 0.6, 0.8, 0.8],),
        ]),

    ("POASU 4", "EG_regime_stableanoxic_to_unstable4", 
        [
            ("k_O2_U", (0.23, 0.25)), 
            ("corg_burial_fac", 1.0), 
            ("OPAinit", (0.3*3.7e19, 4*3.1e15, 4.2*3.193e18)),

            ("set_interp_forcing", "ocean", "force_corg_bf", [-1e30, 5e6 - 1, 5e6, 26e6, 26e6+1, 35e6], [0.3, 0.3, 0.3, 0.3, 0.3, 0.3],),
        ]),
]

isdir(output_figures_dir) || mkdir(output_figures_dir)

P_O_A_S_U_columns_EG = Dict()

for (expt_id, fileroot, vector_pars) in expts_table

# (expt_id, fileroot, vector_pars) = expts_table[1] # yamls, modelname, modelpars,

    ooeoae_expts(
        model, vector_pars
    )

    tstart = 0
    tend =  35e6
    tspan = (tstart, tend) # yr # tspan=(-1000e6, 0)

    #########################################################
    # Initialize
    #########################################################

    initial_state, modeldata = PALEOmodel.initialize!(model)

    ################################
    # get time_series
    ####################################
    # function find_time_series is too slow, back to use the original one
    (paleorun, A_ts, O_ts, P_ts) = SolverFunctionsOOEOAE2.find_time_series(model, initial_state, modeldata; has_A=true, tspan)

    # # plot time series
    filename = "EG_C_U_S.xlsx"
    df_C = DataFrames.DataFrame(XLSX.readtable(filename,"EG-C"))
    df_U = DataFrames.DataFrame(XLSX.readtable(filename,"EG-U"))
    df_S = DataFrames.DataFrame(XLSX.readtable(filename,"EG-S"))

    C_delta = PB.get_data(paleorun.output, "ocean.DIC_delta")
    U_delta = PB.get_delta.(PB.get_data(paleorun.output, "ocean.U"))
    S_delta = PB.get_data(paleorun.output, "ocean.S_delta")
    P_norm = PB.get_data(paleorun.output, "ocean.P_norm")
    U_norm = PB.get_data(paleorun.output, "ocean.U_norm")
    S_norm = PB.get_data(paleorun.output, "ocean.S_norm")
    O_norm = PB.get_data(paleorun.output, "atmocean.O_norm")
    A_norm = PB.get_data(paleorun.output, "atmocean.A_norm")
    pO2atm =  PB.get_data(paleorun.output, "atm.pO2atm")
    pCO2PAL =  PB.get_data(paleorun.output, "atm.pCO2PAL")
    TEMP = PB.get_data(paleorun.output, "global.TEMP")
    corg_bf = PB.get_data(paleorun.output, "ocean.corg_bf")
    ANOX = PB.get_data(paleorun.output, "oceanfloor.anoxia_burial_frac")
    ts = (PB.get_data(paleorun.output, "global.tforce") .- 600e6)/-1e6
    FCorgb = PB.get_total.(PB.get_data(paleorun.output, "fluxOceanBurial.flux_total_Corg"))

    plot_corg_bf = plot(ts, corg_bf, ylabel="Corg burial \n factor", xlabel="", color=:blue, label=false, xlims=(565, 600), xflip=true, ylims=(0, 1.5), yflip=false)
    
    plot_norms = plot(ts, P_norm, label="P_norm"); plot!(plot_norms, ts, O_norm, label="O_norm"); plot!(plot_norms, ts, A_norm, label="A_norm"); plot!(plot_norms, ts, S_norm, label="S_norm");
        plot!(plot_norms, ts, U_norm, ylabel="Reservoirs \n (normalized)", xlabel="", xlims=(565, 600), xflip=true, ylims=(0.0, 6.0), label="U_norm")

    plot()
    # plot(ts, pO2atm, xlabel="", ylabel="pO2atm \n (atm)", color=:blue, yguidefontcolor=:blue, label=false, xlims=(565, 600), xflip=true, ylims=(0, 0.5));
    #     plot_CO2_O2 = plot!(twinx(), ts, pCO2PAL, ylabel="pCO2PAL", color=:red, yguidefontcolor=:red, label=false, xlims=(565, 600), xflip=true, ylims=(10, 24))
    plot_pO2 = plot(ts, pO2atm, xlabel="", ylabel="pO2atm \n (atm)", color=:blue, label=false, xlims=(565, 600), xflip=true, ylims=(0, 0.2))
    plot_pCO2 = plot(ts, pCO2PAL, xlabel="", ylabel="pCO2PAL", color=:blue, label=false, xlims=(565, 600), xflip=true, ylims=ylims=(10, 24))

    plot_TEMP = plot(ts, TEMP,  ylabel="TEMP (k)", xlabel="", color=:blue, label=false, xlims=(565, 600), xflip=true, ylims=(290.0, 295.0), yflip=false)

    plot_EG_C = scatter(df_C.Age, df_C.d13Ccarb, color=:gray, alpha=0.006, ylabel="δ¹³Ccarb", label=false,);
        plot!(plot_EG_C, ts, C_delta, xlabel="tmodel (Ma)", ylabel="δ¹³Ccarb", label=false, color=:blue, xlims=(565, 600), xflip=true, ylims=(-8, 8))

    plot_EG_U = scatter(df_U.Age, df_U.d238U, color=:gray, alpha=0.07, ylabel = "δ²³⁸Ucarb", label=false,)
        plot!(plot_EG_U, ts, U_delta, xlabel="tmodel (Ma)", ylabel="δ²³⁸Ucarb", label=false, color=:blue, ylims=(-1.0, 0.6), xlims=(565, 600), xflip=true,)
    
    plot_EG_S = scatter(df_S.Age, df_S.d34Scas, color=:gray, alpha=0.05, ylabel = "δ³⁴Scas", label=false,)
        plot!(plot_EG_S, ts, S_delta, xlabel="tmodel (Ma)", ylabel="δ³⁴Scas", xlims=(565, 600), xflip=true,)
    
    plot_ANOX = plot(ts, ANOX, xlabel="", ylabel="Degree of anoxia", label=false, color=:blue, xlims=(565, 600), xflip=true,)
    
    
    l = @layout[
            grid(4,2);
    ]

    p_sum_timeseries = Plots.plot(
                                plot_corg_bf, plot_norms, 
                                plot_pCO2,  plot_pO2, 
                                plot_TEMP, plot_ANOX,
                                plot_EG_C, plot_EG_U,
                        layout = l, size=(600, 600)) # left_margin = 5Plots.mm, right_margin = 5Plots.mm,

    savefig(p_sum_timeseries, joinpath(output_figures_dir, "POASU_"*fileroot*"_time_series_20240601.svg"))  

    ##########################################
    # save everything in 3D for the sum plot!
    P_O_A_S_U_columns_EG[fileroot] = 
        (; expt_id, paleorun)
    ##########################################

end # end of for loop

# linestyles = [:solid, :dash, :dashdot]
# colors = [:red, :grey, :grey]

# function plot_3timeseries(var::String)
#     pp = plot()
#     for (i, fileroot) in enumerate(["EG_regime_stableanoxic_to_unstable1", "EG_regime_stableanoxic_to_unstable2", "EG_regime_stableanoxic_to_unstable3"])
#         (; expt_id, paleorun) = P_O_A_S_U_columns_EG[fileroot]
#         data = PB.get_data(paleorun.output, var)
#         ts = (PB.get_data(paleorun.output, "global.tforce") .- 600e6)/-1e6

#         plot!(pp, ts, data, color=colors[i], linestyle=linestyles[i])
#     end
#     return pp
# end

# plot_EG_C = plot_3timeseries("ocean.DIC_delta")


