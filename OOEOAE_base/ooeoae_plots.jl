import Plots
import MAT
import StatsPlots
import PALEOboxes as PB
import PALEOmodel

using Plots, MAT
# using Plots, MAT, StatsPlots
# using Animations

function plot_OOEOAE_columns(
    paleorun,
    fileroot,
    Pbal_P,
    Pbal_O,
    Obal_P, 
    Obal_O,
    P_ts,
    O_ts,
    plot_A_U::String="no_A_U",
)
    if plot_A_U == "no_A_U"
        gr(size=(1400, 300))
        pager=PALEOmodel.PlotPager((1, 4), (legend_background_color=nothing, ); displayfunc=(plot, nplot)->savefig(plot, "plot_$(fileroot).png")) # save to file instead of default display)

        pager(
            (
                scatter(Pbal_P[1].elements, Pbal_O[1].elements, xlims=[0.0, 5.5], ylims=[0.0, 2.25], label="slow manifold", xlabel="P_norm", ylabel="O (PAL)");
                plot!(Pbal_P[1].elements, Pbal_O[1].elements, xlims=[0.0, 5.5], ylims=[0.0, 2.25], label="slow manifold", xlabel="P_norm", ylabel="O (PAL)");
                plot!(Obal_P[1].elements, Obal_O[1].elements, label="O nullcline");
                plot!(P_ts, O_ts, label="time series", left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
            ),
            plot(ylabel="Reservoirs", paleorun.output, ["atmocean.O_norm", "ocean.P_norm", "sedcrust.C_norm", "sedcrust.G_norm"]),

            plot(paleorun.output, ylabel="P fluxes", ["fluxRtoOcean.flux_P", "fluxOceanBurial.flux_total_P"]), 
            plot(paleorun.output, ylabel="O fluxes", ["fluxOceanBurial.flux_total_Corg", "fluxAtoLand.flux_O2", "fluxSedCrusttoAOcean.flux_Redox"]), 
        )
    elseif plot_A_U == "A_U"
        gr(size=(1400, 600)) # plotlyjs(size=(1600, 900))
        pager=PALEOmodel.PlotPager((2, 4), (legend_background_color=nothing, ); displayfunc=(plot, nplot)->savefig(plot, "plot_$(fileroot).png")) # save to file instead of default display)

        pager(
            (
                scatter(Pbal_P[1].elements, Pbal_O[1].elements, xlims=[0.0, 5.5], ylims=[0.0, 2.25], label="slow manifold", xlabel="P_norm", ylabel="O (PAL)");
                plot!(Pbal_P[1].elements, Pbal_O[1].elements, xlims=[0.0, 5.5], ylims=[0.0, 2.25], label="slow manifold", xlabel="P_norm", ylabel="O (PAL)");
                plot!(Obal_P[1].elements, Obal_O[1].elements, label="O nullcline");
                plot!(P_ts, O_ts, label="time series", left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
            ),
            plot(ylabel="Reservoirs", paleorun.output, ["atmocean.O_norm", "ocean.P_norm", "sedcrust.C_norm", "sedcrust.G_norm", "atmocean.A_norm", "ocean.U_norm"]),
            plot(paleorun.output, ylabel="P fluxes", ["fluxRtoOcean.flux_P", "fluxOceanBurial.flux_total_P"]), 
            plot(paleorun.output, ylabel="O fluxes", ["fluxOceanBurial.flux_total_Corg", "fluxAtoLand.flux_O2", "fluxSedCrusttoAOcean.flux_Redox"]), 

            plot(ylabel="TEMP (K)", paleorun.output, ["global.TEMP"]),
            plot(title="Carbon isotopes",  paleorun.output, ["atmocean.A_delta", "ocean.DIC_delta", "atm.CO2_delta"], ylabel="delta 13C (per mil)"), # , xlims=(-1e6, 10e6), "ocean.mccb_delta", "sedcrust.C_delta" ; extrakwargs...),
            plot(title="Uranium isotopes",  paleorun.output, ["fluxRtoOcean.flux_U.v_delta", "ocean.U.v_delta"], ylabel="d238U/235U"),
            plot(title="anoxia_burial_frac", paleorun.output, ["oceanfloor.anoxia_burial_frac"])
        )
    end

    pager(:newpage)

end


"""
    plot_anim_P_O

    For the animation plots
"""
function plot_anim_P_O(
    Pbal_P,
    Pbal_O,
    Obal_P, 
    Obal_O,
    P_ts,
    O_ts,
    plot_name::String="default",
)
    num_O_ts = length(O_ts)
    gr(size=(500, 500))
    pp = Plots.plot(Pbal_P, Pbal_O, xlims=[0.0, 2.5], ylims=[0.0, 1.25], label="slow manifold", xlabel="P_norm", ylabel="O (PAL)")
    pp = Plots.plot!(Obal_P, Obal_O, label="O nullcline")

    function get_Ots(i)
        if i in 1:num_O_ts
            return O_ts[i]
        else
            return i
        end
    end

    # pp = scatter!([get_Ots], zeros(0), lw=3, xlims=[0.0, 2.5], ylims=[0.0, 1.25])
    # anim = @animate for i in eachindex(P_ts)
    #     push!(pp, P_ts[i], Float64[get_Ots(i)])
    # end

    anim = @animate for i in eachindex(P_ts)
        Plots.scatter!(pp, [P_ts[i]], [get_Ots(i)], xlims=[0.0, 5.5], ylims=[0.0, 2.25], color=:green, label=false)
    end

    if plot_name=="default"
        gif(anim, fps=20)
    else
        gif(anim, plot_name, fps=20)
    end
end




