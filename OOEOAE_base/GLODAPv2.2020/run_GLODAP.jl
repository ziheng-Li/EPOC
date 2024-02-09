import Statistics
import DataFrames
import MAT
import StatsPlots
import NCDatasets

using DataFrames, StatsPlots, NCDatasets

include("base_GLODAP.jl")

# filenameroots = ["Arctic_Ocean", "Atlantic_Ocean", "Indian_Ocean", "Pacific_Ocean", "Merged_master"]

# df = read_GLODAP(filenameroots[2])

df_shelf = Data_filter(["Arctic_Ocean" "Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean"], 6000.0, "all", "oxygen", 10)

histogram2d(df_shelf.G2phosphate, df_shelf.G2depth, bins=(40, 40), yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,6000], color=:plasma)

# plot raw data
# plot_TEMP_O2_P_dC13(["Arctic_Ocean"])

# three options: polar ocean, inner ocean, epicontinental, following GeoClim
# plot_TEMP_O2_P_dC13(["Atlantic_Ocean"  "Indian_Ocean"  "Pacific_Ocean"], "highlat", "oxygen")
# plot_TEMP_O2_P_dC13(["Atlantic_Ocean"  "Indian_Ocean"  "Pacific_Ocean"], "lowlat", "oxygen")
# plot_TEMP_O2_P_dC13(["Atlantic_Ocean"  "Indian_Ocean"  "Pacific_Ocean"], 200.0, "all", "oxygen")





