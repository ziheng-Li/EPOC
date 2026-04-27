using Logging
import DataFrames
import Interpolations
using Plots
import XLSX

import PALEOboxes as PB
import PALEOmodel
# import PALEOreactions
import PALEOcopse

include("../../src/EPOC_reactions.jl")
include("ooeoae_expts.jl")
include("copse_expts.jl")

######################################
# # key parameters that control plot appearance
# include("expt_plot3D_EG.jl")
P_lims = (0.0, 5.0)
O_lims = (0.0, 2.0)
A_lims = (2, 5)

dropbox_output_dir = "/Users/liziheng/Dropbox/SPICE"
output_figures_dir = joinpath(dropbox_output_dir, "Figs")
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
        ("CPsea", 115.4, 115.4*3), # base setup values for Bergman_COPSE
        ("P_weathering", 3.9e10, 1.0, 0.0, 0.0),
        ("CO2pulse", [0, 5e7, 7.0e7, 1e8], [0, 0, 0, 0]), 
        ("COPSE_locb", 1000.0, 0.0),
        ("k_anox", 100.0),
    ]
)

# Experiments as listed in Table2 in README
expts_table = [
    # ("POASU 1", "EG_regime_stableanoxic_to_unstable1", (490, 497),
    #     [
    #         ("k_O2_U", (0.23, 0.25)), 
    #         ("corg_burial_fac", 1.0), 
    #         ("OPAinit", (0.3*3.7e19, 4*3.1e15, 4.2*3.193e18)),

    #         ("set_interp_forcing", "ocean", "force_corg_bf", [-1e30, 5e6 - 1, 5e6, 26e6, 26e6+1, 35e6], [0.8, 0.8, 0.8, 0.8, 0.8, 0.8],),
    #     ]),
    
    ("POASU 1", "EG_regime_stableanoxic_to_unstable1", (490, 540),
        [
            ("k_O2_U", (0.23, 0.25)), 
            ("corg_burial_fac", 1.0), 
            ("OPAinit", (0.3*3.7e19, 4*3.1e15, 4.2*3.193e18)),

            ("set_interp_forcing", "ocean", "force_corg_bf", [-1e30, 20e6, 20e6+1, 55e6], [0.8, 0.8, 1.0, 1.0],),
        ]),
    
    # ("POASU 2", "EG_regime_stableanoxic_to_unstable2", 
    #     [
    #         ("k_O2_U", (0.23, 0.25)), 
    #         ("corg_burial_fac", 1.0), 
    #         ("OPAinit", (0.3*3.7e19, 4*3.1e15, 4.2*3.193e18)),

    #         ("set_interp_forcing", "ocean", "force_corg_bf", [-1e30, 5e6 - 1, 5e6, 26e6, 26e6+1, 35e6], [0.3, 0.3, 0.8, 0.8, 1.1, 1.1],),
    #     ]),

    # ("POASU 3", "EG_regime_stableanoxic_to_unstable3", 
    #     [
    #         ("k_O2_U", (0.23, 0.25)), 
    #         ("corg_burial_fac", 1.0), 
    #         ("OPAinit", (0.3*3.7e19, 4*3.1e15, 4.2*3.193e18)),

    #         ("set_interp_forcing", "ocean", "force_corg_bf", [-1e30, 5e6 - 1, 5e6, 26e6, 26e6+1, 35e6], [0.3, 0.3, 0.6, 0.6, 0.8, 0.8],),
    #     ]),

    # ("POASU 4", "EG_regime_stableanoxic_to_unstable4", 
    #     [
    #         ("k_O2_U", (0.23, 0.25)), 
    #         ("corg_burial_fac", 1.0), 
    #         ("OPAinit", (0.3*3.7e19, 4*3.1e15, 4.2*3.193e18)),

    #         ("set_interp_forcing", "ocean", "force_corg_bf", [-1e30, 5e6 - 1, 5e6, 26e6, 26e6+1, 35e6], [0.3, 0.3, 0.3, 0.3, 0.3, 0.3],),
    #     ]),
]

isdir(output_figures_dir) || mkdir(output_figures_dir)

P_O_A_S_U_columns_SPICE = Dict()
(plot_SPICE_C, plot_SPICE_U, plot_SPICE_S, df_C, df_U, df_S) = plot_C_U_S("/Users/liziheng/Dropbox/SPICE/Spice-20250522.xlsx")


# for (expt_id, fileroot, xlims, vector_pars) in expts_table

(expt_id, fileroot, xlims, vector_pars) = expts_table[1] # yamls, modelname, modelpars,

    ooeoae_expts(
        model, vector_pars
    )

    tstart = 0
    tend =  55e6
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
    filename = "/Users/liziheng/Dropbox/SPICE/Spice-20250522.xlsx"
    df_C = DataFrames.DataFrame(XLSX.readtable(filename,"C"))
    df_U = DataFrames.DataFrame(XLSX.readtable(filename,"U"))
    df_S = DataFrames.DataFrame(XLSX.readtable(filename,"S"))

    C_delta = PB.get_data(paleorun.output, "ocean.DIC_delta")
    U_delta = PB.get_delta.(PB.get_data(paleorun.output, "ocean.U"))
    S_delta = PB.get_data(paleorun.output, "ocean.S_delta")
    P_norm = PB.get_data(paleorun.output, "ocean.P_norm")
    U_norm = PB.get_data(paleorun.output, "ocean.U_norm")
    S_norm = PB.get_data(paleorun.output, "ocean.S_norm")
    O_norm = PB.get_data(paleorun.output, "atmocean.O_norm")
    A_norm = PB.get_data(paleorun.output, "atmocean.A_norm")
    pO2PAL =  PB.get_data(paleorun.output, "atm.pO2PAL")
    pCO2PAL =  PB.get_data(paleorun.output, "atm.pCO2PAL")
    TEMP = PB.get_data(paleorun.output, "global.TEMP")
    corg_bf = PB.get_data(paleorun.output, "ocean.corg_bf")
    ANOX = PB.get_data(paleorun.output, "oceanfloor.anoxia_burial_frac")
    ts = (PB.get_data(paleorun.output, "global.tforce") .- 540.5e6)/-1e6
    FCorgb = PB.get_total.(PB.get_data(paleorun.output, "fluxOceanBurial.flux_total_Corg"))

    ######## for the plots of source and sinks ########
    Fmocb = PB.get_total.(PB.get_data(paleorun.output, "fluxOceanBurial.flux_total_Corg"))
    Flocb = PB.get_total.(PB.get_data(paleorun.output, "land.locb"))
    Fmccb_percell = PB.get_total.(PB.get_data(paleorun.output, "oceanfloor.mccb"))
    Fmccb=[]
    for i = 1:length(Fmccb_percell)
        push!(Fmccb, sum(Fmccb_percell[i]))
    end
    Foxidw = PB.get_data(paleorun.output, "land.oxidw")
    Fcarbw = PB.get_data(paleorun.output, "land.carbw")
    Focdeg = PB.get_data(paleorun.output, "sedcrust.ocdeg")
    Fccdeg = PB.get_data(paleorun.output, "sedcrust.ccdeg")
    Fpsea = PB.get_data(paleorun.output, "land.psea")
    Fmpb = PB.get_data(paleorun.output, "fluxOceanBurial.flux_total_P")
    Fpyrw = PB.get_total.(PB.get_data(paleorun.output, "land.S_pyrite_weath"))
    Fpyrb = PB.get_total.(PB.get_data(paleorun.output, "oceanfloor.S_pyrite_total"))
    Fuanoxic = PB.get_total.(PB.get_data(paleorun.output, "oceanfloor.U_anoxic"))
    Fuother = PB.get_total.(PB.get_data(paleorun.output, "oceanfloor.U_other"))
    ################################################################

    plot_corg_bf = plot(ts, corg_bf, ylabel="Corgbf", xlabel="", color=:blue, label=false, xlims=xlims, xflip=true, ylims=(0, 1.5), yflip=false)
    
    plot_norms = plot(ts, P_norm, label="P_norm"); plot!(plot_norms, ts, O_norm, label="O_norm"); plot!(plot_norms, ts, A_norm, label="A_norm"); plot!(plot_norms, ts, S_norm, label="S_norm");
        plot!(plot_norms, ts, U_norm, ylabel="Reservoirs \n (normalized)", xlabel="tmodel (Ma)", xlims=xlims, xflip=true, ylims=(0.0, 4.5), label="U_norm", legend=:topright)

    plot()
    # plot(ts, pO2PAL, xlabel="", ylabel="pO2PAL \n (atm)", color=:blue, yguidefontcolor=:blue, label=false, xlims=xlims, xflip=true, ylims=(0, 0.5));
    #     plot_CO2_O2 = plot!(twinx(), ts, pCO2PAL, ylabel="pCO2PAL", color=:red, yguidefontcolor=:red, label=false, xlims=xlims, xflip=true, ylims=(10, 24))
    plot_pO2 = plot(ts, pO2PAL, xlabel="", ylabel="pO2 (PAL)", color=:blue, label=false, xlims=xlims, xflip=true, ylims=(0.2, 0.6))
    plot_pCO2 = plot(ts, pCO2PAL, xlabel="", ylabel="pCO2PAL", color=:blue, label=false, xlims=xlims, xflip=true, ylims=ylims=(10, 24))

    plot_TEMP = plot(ts, TEMP,  ylabel="TEMP (k)", xlabel="", color=:blue, label=false, xlims=xlims, xflip=true, ylims=(290.0, 295.0), yflip=false)

    plot_EG_C = scatter(df_C.Age*-1, df_C.d13C, color=:gray, alpha=0.02, ylabel="δ¹³Ccarb", label=false,);
        plot!(plot_EG_C, ts, C_delta, ylabel="δ¹³Ccarb", label=false, color=:blue, xlims=xlims, xflip=true, ylims=(-3, 6))

    plot_EG_U = scatter(df_U.Age*-1, df_U.d238U, color=:gray, alpha=0.08, ylabel = "δ²³⁸Ucarb", label=false,)
        plot!(plot_EG_U, ts, U_delta, ylabel="δ²³⁸Ucarb", label=false, color=:blue, ylims=(-1.0, -0.2), xlims=xlims, xflip=true,)
    
    plot_EG_S = scatter(df_S.Age*-1, df_S.d34SCAS, color=:gray, alpha=0.08, ylabel = "δ³⁴Scas", label=false,)
        plot!(plot_EG_S, ts, S_delta, xlabel="tmodel (Ma)", ylabel="δ³⁴Scas", label=false, color=:blue, xlims=xlims, xflip=true, ylims=(20, 70))
    
    plot_ANOX = plot(ts, ANOX, xlabel="tmodel (Ma)", ylabel="ANOX", label=false, color=:blue, xlims=xlims, xflip=true,)

    #######################################
    plot_A_ss = plot(ts, (Foxidw+Fcarbw+Focdeg+Fccdeg)/1e13, label="A source", ylabel="A sourve & sink \n (x 1e13 mol/yr)", xlims=(490, 497), ylims=(1.5, 3.0), xflip=true)
        plot!(plot_A_ss, ts, (Fmocb+Fmccb)/1e13, label="A sink")

    plot_P_ss = plot(ts, Fpsea/1e10, label="P source", ylabel="P sourve & sink \n (x 1e10 mol/yr)", xlims=(490, 497), ylims=(1.6, 7.0), xflip=true)
        plot!(plot_P_ss, ts, Fmpb/1e10, label="P sink", xlabel="tmodel (Ma)")
    
    plot_O_ss = plot(ts, (Fmocb+Flocb+2*Fpyrb)/1e13, label="O source", ylabel="O sourve & sink \n (x 1e13 mol/yr)", xlims=(490, 497), ylims=(0.4, 1.6), xflip=true)
        plot!(plot_O_ss, ts, (Foxidw+Focdeg+2*Fpyrw)/1e13, label="O sink")
    
    plot_S_ss = plot(ts, Fpyrw, label="S source", ylabel="S sourve & sink \n (x 1e13 mol/yr)", xlims=(490, 497),  xflip=true) # ylims=(0.4, 1.6),
        plot!(plot_S_ss, ts, Fpyrb, label="S sink")
    
    plot_pO2_zoomin = plot(ts, pO2PAL, xlabel="", ylabel=false, color=:red, label="pO₂", xlims=(490, 497), xflip=true, ylims=(-0.1, 1.1))
        plot!(plot_pO2_zoomin, ts, ANOX, ylabel="fANOX & \n pO₂ (PAL)", label="fANOX", color=:blue, xlims=(490, 497), xflip=true,)
    
    plot_norms_zoomin = plot(ts, P_norm, label="P_norm"); plot!(plot_norms_zoomin, ts, O_norm, label="O_norm"); plot!(plot_norms_zoomin, ts, A_norm, label="A_norm"); plot!(plot_norms_zoomin, ts, S_norm, label="S_norm");
        plot!(plot_norms_zoomin, ts, U_norm, ylabel="Reservoirs \n (normalized)", xflip=true, ylims=(0.01, 10), xlims=(490, 497), label="U_norm", yaxis=:log, legend=:topright)

    plot_corg_bf_zoomin = plot(ts, corg_bf, ylabel="Corgbf", xlabel="", color=:blue, label=false, xflip=true, ylims=(0, 1.5), xlims=(490, 497), yflip=false)

    plot_EG_C_zoomin = scatter(df_C.Age*-1, df_C.d13C, color=:gray, alpha=0.06, ylabel="δ¹³Ccarb", label=false,);
        plot!(plot_EG_C_zoomin, ts, C_delta, ylabel="δ¹³Ccarb (‰)", label=false, color=:blue, xlims=(490, 497), xflip=true, ylims=(-3, 6))

    plot_EG_U_zoomin = scatter(df_U.Age*-1, df_U.d238U, color=:gray, alpha=0.1, ylabel = "δ²³⁸Ucarb", label=false,)
        plot!(plot_EG_U_zoomin, ts, U_delta, ylabel="δ²³⁸Ucarb (‰)", label=false, color=:blue, ylims=(-1.0, -0.2), xlims=(490, 497), xflip=true,)
    
    plot_EG_S_zoomin = scatter(df_S.Age*-1, df_S.d34SCAS, color=:gray, alpha=0.1, ylabel = "δ³⁴SCAS", label=false,)
        plot!(plot_EG_S_zoomin, ts, S_delta, xlabel="tmodel (Ma)", ylabel="δ³⁴SCAS (‰)", label=false, color=:blue, xlims=(490, 497), xflip=true, ylims=(20, 70))
    

    l = @layout[
        [a{0.48w} grid(2,1)]
        
        grid(2,2)
    ]

    title=["a","b","c","d","e","f","g","h","i","J","K", "L"]

    p_sum_timeseries = Plots.plot(
        plot_O_ss, plot_pO2_zoomin, plot_EG_C_zoomin,
                                plot_A_ss, plot_EG_U_zoomin,
                                plot_P_ss, plot_EG_S_zoomin,
                        layout = l, size=(600, 700), right_margin = 2Plots.mm,
                        title = ["$(title[i])" for j in 1:1, i in 1:8], titleloc = :left, titlefont = font(10)
     ) # left_margin = 2Plots.mm,

    display(p_sum_timeseries)
    savefig(p_sum_timeseries, joinpath(output_figures_dir, "POASU_"*fileroot*"_time_series_20250810.svg"))  

    ##########################################
    # save everything in 3D for the sum plot!
    P_O_A_S_U_columns_SPICE[fileroot] = 
        (; expt_id, paleorun)
    ##########################################

# end # end of for loop

# linestyles = [:solid, :dash, :dashdot]
# colors = [:red, :grey, :grey]

# function plot_3timeseries(var::String)
#     pp = plot()
#     for (i, fileroot) in enumerate(["EG_regime_stableanoxic_to_unstable1", "EG_regime_stableanoxic_to_unstable2", "EG_regime_stableanoxic_to_unstable3"])
#         (; expt_id, paleorun) = P_O_A_S_U_columns_SPICE[fileroot]
#         data = PB.get_data(paleorun.output, var)
#         ts = (PB.get_data(paleorun.output, "global.tforce") .- 600e6)/-1e6

#         plot!(pp, ts, data, color=colors[i], linestyle=linestyles[i])
#     end
#     return pp
# end

# plot_EG_C = plot_3timeseries("ocean.DIC_delta")


