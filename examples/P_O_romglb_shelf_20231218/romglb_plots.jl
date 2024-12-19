import Plots
import MAT
import StatsPlots
import PALEOboxes as PB
import PALEOmodel

using Plots, MAT
# using Plots, MAT, StatsPlots
# using Animations


function plot_ocean_tracers(
    output;
    columns=[:hlat, :gyre, :upw],
    colskip=0, # number of blank panels to add (use to format pages)
    tracers=[
        "insol", "DIC_conc", "TAlk_conc", "temp", "pHtot", "O2_conc", "SO4_conc", "H2S_conc", "CH4_conc","P_conc", "Ca_conc",
        "H2S_delta", "SO4_delta", "CH4_delta", "OmegaAR"],
    tcol=[-Inf, 10.0, 100.0, 1000.0, Inf], # model time for column plots
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)
    for tr in tracers
        for col in columns
            pager(plot(title="$tr $col", output, "ocean.$tr", ( tmodel=tcol, column=col);
                    swap_xy=true, labelattribute=:filter_records, plotargs...))
        end
        for i in 1:colskip
            pager(:skip)
        end
    end

    return nothing
end

function plot_Corg_burial_alongdepth(
    output,
)

    corg_burial = PB.get_data(output, "fluxOceanBurial.flux_Corg")
    ocean_Z_mid = PB.get_data(output, "ocean.zmid")

    p1 = plot(corg_burial[end][1:13],ocean_Z_mid[end][1:13], label="hlat")
    plot!(corg_burial[end][14:46],ocean_Z_mid[end][14:46], label="gyre")
    plot!(corg_burial[end][47:82],ocean_Z_mid[end][47:82], label="shelves")

    depth=collect(0:50:1000); append!(depth, collect(2000:1000:6000))

    display(p1)
end

"""
    plot_P_vs_burial


"""
function plot_P_vs_burial(grid_P_norm, grid_pO2PAL, burialflux, ylabel="No title")

    pp = plot()
    for i in eachindex(grid_pO2PAL)
        plot!(grid_P_norm, burialflux[i,:], label="pO2PAL=$(grid_pO2PAL[i])", xlabel="P_norm", ylabel=ylabel, left_margin = 5Plots.mm, right_margin = 5Plots.mm)
    end

    return pp
end

"""
    plot_compare_GLODAP

    Make the comparison of Romanielle's ocean and GLODAP's data
"""
function plot_compare_GLODAP(
    output,
    df_hlat,
    df_gyre,
    df_upw,
    df_shelf,
    savefig_name = "NaN",
    )

    gr(size=(1200, 600))

    # df_hlat = Data_filter(["Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 6000.0, "highlat", "oxygen")
    # df_gyre = Data_filter(["Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 6000.0, "lowlat", "oxygen")
    # df_upw = Data_filter(["Arctic_Ocean" "Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 6000.0, "all", "oxygen", 10)
    # df_shelf = Data_filter(["Arctic_Ocean" "Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 200.0, "all", "oxygen")

    O2_conc = PB.get_data(output, "ocean.O2_conc") *1000 # convert to mmmol/m3
    P_conc = PB.get_data(output, "ocean.P_conc") *1000 # convert to mmmol/m3
    ocean_depth = PB.get_data(output, "ocean.zmid") *-1 # convert
    Vector_num_cells = PB.get_data(output, "ocean.number_cells")
    num_cells = Int64(Vector_num_cells[1])

    # High latitude plots
    p1 = histogram2d(df_hlat.G2oxygen, df_hlat.G2depth, bins=(40, 40), title="Hlat", yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,6000], color=:YlGnBu_9)
    Plots.plot!(p1, O2_conc[end][1:13], ocean_depth[end][1:13], label="Model", color=:red, linewidth=2, marker=:circle)
    p2 = histogram2d(df_hlat.G2phosphate, df_hlat.G2depth, bins=(40, 40), title="Hlat", yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,6000], color=:YlOrRd_9)
    Plots.plot!(p2, P_conc[end][1:13], ocean_depth[end][1:13], label="Model", color=:red, linewidth=2, marker=:rtriangle)

    # low latitude plots
    p3 = histogram2d(df_gyre.G2oxygen, df_gyre.G2depth, bins=(40, 40), title="Gyre", yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,6000], color=:YlGnBu_9)
    Plots.plot!(p3, O2_conc[end][14:46], ocean_depth[end][14:46], label=false, color=:red, linewidth=2, marker=:circle)
    p4 = histogram2d(df_gyre.G2phosphate, df_gyre.G2depth, bins=(40, 40), title="Gyre", yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,6000], color=:YlOrRd_9)
    Plots.plot!(p4, P_conc[end][14:46], ocean_depth[end][14:46], label=false, color=:red, linewidth=2, marker=:rtriangle)

    # upw column
    if num_cells > 46
        p5 = histogram2d(df_upw.G2oxygen, df_upw.G2depth, bins=(40, 40), title="Upw", yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,6000], color=:YlGnBu_9)
        Plots.plot!(p5, O2_conc[end][47:79], ocean_depth[end][47:79], label=false, color=:red, linewidth=2, marker=:circle)
        p6 = histogram2d(df_upw.G2phosphate, df_upw.G2depth, bins=(40, 40), title="Upw", yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,6000], color=:YlOrRd_9)
        Plots.plot!(p6, P_conc[end][47:79], ocean_depth[end][47:79], label=false, color=:red, linewidth=2, marker=:rtriangle)
    else
        p5 = histogram2d(df_upw.G2oxygen, df_upw.G2depth, bins=(40, 40), title="Upw", yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,6000], color=:YlGnBu_9)
        p6 = histogram2d(df_upw.G2phosphate, df_upw.G2depth, bins=(40, 40), title="Upw", yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,6000], color=:YlOrRd_9)
    end

    if num_cells > 79
        # shelves
        p7 = histogram2d(df_shelf.G2oxygen, df_shelf.G2depth, bins=(40, 40), title="Shelf", yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,200], color=:YlGnBu_9)
        Plots.plot!(p7, O2_conc[end][80:num_cells], ocean_depth[end][80:num_cells], label=false, color=:red, linewidth=2, marker=:circle)
        p8 = histogram2d(df_shelf.G2phosphate, df_shelf.G2depth, bins=(40, 40), title="Shelf", yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,200], color=:YlOrRd_9)
        Plots.plot!(p8, P_conc[end][80:num_cells], ocean_depth[end][80:num_cells], label=false, color=:red, linewidth=2, marker=:rtriangle)
    else
        p7 = histogram2d(df_shelf.G2oxygen, df_shelf.G2depth, bins=(40, 40), title="Shelf", yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,200], color=:YlGnBu_9)
        p8 = histogram2d(df_shelf.G2phosphate, df_shelf.G2depth, bins=(40, 40), title="Shelf", yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,200], color=:YlOrRd_9)
    end

    l = @layout[a b c d; e f g h]
    Plots.plot(p7, p5, p3, p1, p8, p6, p4, p2, layout = l, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)

    if savefig_name != "NaN"
        savefig(savefig_name)
    else
        # display(plot(p1, p2, p3, p4, p5, p6, layout = l))
        display(Plots.plot(p7, p5, p3, p1, p8, p6, p4, p2, layout = l, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm))
    end

end

"""
    plot_compare_GLODAP_shelf_only

    Make the plots of O and P of GLODAP, input a vector of different subset of data.
    Make the comparison of model and GLODAP
"""
function plot_compare_GLODAP_shelf_only(
    output,
    vector_database=[],
    )

    gr(size=(600, 800))

    O2_conc = PB.get_data(output, "ocean.O2_conc") *1000 # convert to mmmol/m3
    P_conc = PB.get_data(output, "ocean.P_conc") *1000 # convert to mmmol/m3
    ocean_depth = PB.get_data(paleorun.output, "ocean.zmid") *-1 # convert
    Vector_num_cells = PB.get_data(paleorun.output, "ocean.number_cells")
    num_cells = Int64(Vector_num_cells[1])

    if length(vector_database) != 4
        @error "The function plot_P_O_GLODAP_shelf needs four subsets of data represent [tropical, subtropical, subpolar, polar] shelves"
    end

    if (num_cells - 46) / 2 != 18
        @error "The function plot_P_O_GLODAP_shelf needs 18 shelves, only have $((num_cells - 46) / 2) right now"
    end

    # Tropical shelf first three shelvess
    p1 = histogram2d(vector_database[1].G2oxygen, vector_database[1].G2depth, bins=(40, 40), title="Tropical shelf", yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,200], color=:YlGnBu_9)
    plot!(p1, O2_conc[end][47:52], ocean_depth[end][47:52], label= false, color=:red, linewidth=2, marker=:circle)
    p2 = histogram2d(vector_database[1].G2phosphate, vector_database[1].G2depth, bins=(40, 40), title="Tropical shelf", yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,200], color=:YlOrRd_9)
    plot!(p2, P_conc[end][47:52], ocean_depth[end][47:52], label= false, color=:red, linewidth=2, marker=:rtriangle)

    # sub Tropical shelf
    p3 = histogram2d(vector_database[2].G2oxygen, vector_database[2].G2depth, bins=(40, 40), title="subTropical shelf", yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,200], color=:YlGnBu_9)
    plot!(p3, O2_conc[end][53:62], ocean_depth[end][53:62], label= false, color=:red, linewidth=2, marker=:circle)
    p4 = histogram2d(vector_database[2].G2phosphate, vector_database[2].G2depth, bins=(40, 40), title="subTropical shelf", yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,200], color=:YlOrRd_9)
    plot!(p4, P_conc[end][53:62], ocean_depth[end][53:62], label= false, color=:red, linewidth=2, marker=:rtriangle)

    # sub polar shelf
    p5 = histogram2d(vector_database[3].G2oxygen, vector_database[3].G2depth, bins=(40, 40), title="subPolar shelf", yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,200], color=:YlGnBu_9)
    plot!(p5, O2_conc[end][63:74], ocean_depth[end][63:74], label= false, color=:red, linewidth=2, marker=:circle)
    p6 = histogram2d(vector_database[3].G2phosphate, vector_database[3].G2depth, bins=(40, 40), title="subPolar shelf", yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,200], color=:YlOrRd_9)
    plot!(p6, P_conc[end][63:74], ocean_depth[end][63:74], label= false, color=:red, linewidth=2, marker=:rtriangle)

    # polar
    p7 = histogram2d(vector_database[4].G2oxygen, vector_database[4].G2depth, bins=(40, 40), title="Polar shelf", yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,200], color=:YlGnBu_9)
    plot!(p7, O2_conc[end][75:82], ocean_depth[end][75:82], label= false, color=:red, linewidth=2, marker=:circle)
    p8 = histogram2d(vector_database[4].G2phosphate, vector_database[4].G2depth, bins=(40, 40), title="Polar shelf", yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,200], color=:YlOrRd_9)
    plot!(p8, P_conc[end][75:82], ocean_depth[end][75:82], label= false, color=:red, linewidth=2, marker=:rtriangle)
    

    l = @layout[a b; c d; e f; g h]

    display(plot(p1, p2, p3, p4, p5, p6, p7, p8, layout = l))

end

"""
plot_compare_vs_depth

Make the comparison of the physical pars between model and Romglb's file 
"""
function plot_compare_vs_depth(
    output,
    savefig_name = "NaN",
)
    gr(size=(1000, 600))

    # get data from the model
    ocean_area = PB.get_data(output, "ocean.Abox")
    ocean_zupper = PB.get_data(output, "ocean.zupper")
    ocean_zmid = PB.get_data(output, "ocean.zmid")
    ocean_zlower = PB.get_data(output, "ocean.zlower")
    thickness_model = ocean_zupper .- ocean_zlower
    ocean_volume_cal = PB.get_data(output, "ocean.volume")

    # get data from romglb's mat file
    matdir = joinpath(PALEOreactions.srcdir(), "ocean")
    matfilename = joinpath(matdir, "romaniello_global79.mat")
    rom_data = MAT.matread(matfilename)
    V = vec(rom_data["V"]) 
    box_thickness = vec(rom_data["box_thickness"])
    Hyps = rom_data["Hyps"]; Hyps_79 = vcat(Hyps[1:13], Hyps[15:47], Hyps[49:81])

    # make plots
    p1 = plot_vs_depth(output, ocean_area[1], "Ocean area (m2)", "Depth (m)")
    p2 = plot_vs_depth(output, Hyps_79, "Hpys .mat (m2)", "Depth (m)")
    p3 = plot_vs_depth(output, (ocean_area[1] .- Hyps_79)./(Hyps_79), "Diff, Cal-Hpys", "Depth (m)", (-0.5, 1.5))

    p4 = plot_vs_depth(output, thickness_model[1], "Box thickness (m)", "Depth (m)")
    p5 = plot_vs_depth(output, box_thickness, "box_thickness .mat (m)", "Depth (m)")
    p6 = plot_vs_depth(output, (thickness_model[1] .- box_thickness)./(box_thickness), "Diff, Cal-BT", "Depth (m)", (-0.5, 1.5))

    p7 = plot_vs_depth(output, ocean_volume_cal[1], "Ocean volume calculated (m3)", "Depth (m)")
    p8 = plot_vs_depth(output, V/1000, "V .mat (m3)", "Depth (m)")
    p9 = plot_vs_depth(output, (ocean_volume_cal[1] .- V/1000)./(V/1000), "Diff, Cal-V", "Depth (m)", (-0.5, 1.5))

    l = @layout[a b c; d e f; g h i]
    display(plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout = l))

    if savefig_name != "NaN"
        savefig(savefig_name)
    end
end

"""
plot_vs_depth

Make the cross plot of depth and the input data, 
e.g. the input data could be the bioproduction values per-cell of the last time step of the model output 
"""
function plot_vs_depth(output, data_name, xlabel::String, ylabel::String, xlims=NaN, title="", showdetails=false)
    ocean_zupper = PB.get_data(output, "ocean.zupper")
    ocean_zmid = PB.get_data(output, "ocean.zmid")
    ocean_zlower = PB.get_data(output, "ocean.zlower")

    if typeof(data_name) <: AbstractString
        data = PB.get_data(output, data_name)[end]
    elseif typeof(data_name) <: AbstractArray
        data = data_name
    end

    # # Box indices:
    # #   - hlat      1-13  (surface box is 1)
    # #   - gyre      14-46 
    # #   - upw       47-79
    # #   - shelfX    80-125

    halt_p=""; gyre_p=""; upw_p=""; shelves_p=""; 
    if data_name == "ocean.Abox"
        # if we need the proportion of the surface area, then its different
        # we only need the surface box only!
        if showdetails # then show the proportion of data in each column
            sum_data = data[1] + data[14] + data[47] + sum(data[80:2:end])
            halt_p = conv_to_Str(data[1]/sum_data)
            gyre_p = conv_to_Str(data[14]/sum_data)
            upw_p = conv_to_Str(data[47]/sum_data)
            if length(data) > 79
                shelves_p = conv_to_Str(sum(data[80:2:end])/sum_data)
            end
            @info "The total sum of $(data_name) is $(sum_data)"
        end
    else
        if showdetails # then show the proportion of data in each column
            halt_p = conv_to_Str(sum(data[1:13])/sum(data))
            gyre_p = conv_to_Str(sum(data[14:46])/sum(data))
            upw_p = conv_to_Str(sum(data[47:79])/sum(data))
            if length(data) > 79
                shelves_p = conv_to_Str(sum(data[80:end])/sum(data))
            end
            @info "The total sum of $(data_name) is $(sum(data))"
        end
    end

    pp = plot()
    pp = plot(data[1:13], ocean_zupper[1][1:13], label="hlat"*halt_p)
    pp = plot!(data[14:46], ocean_zupper[1][14:46], label="gyre"*gyre_p)
    pp = plot!(data[47:79], ocean_zupper[1][47:79], label="upw"*upw_p, xlabel=xlabel, xlims=xlims, ylabel=ylabel, title=title, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
    if length(data) > 79
        pp = plot!(data[80:end], ocean_zupper[1][80:end], label="shelves"*shelves_p)
    end
    
    return pp
end

"""
    conv_to_Str

Input e.g. 0.1923567, output " (19.2%)"
"""
function conv_to_Str(Input_float)
    output_Str = string(round(Input_float * 100, digits=1))
    output_Str = " ("*output_Str*"%)"
    return output_Str
end

