import Statistics
import DataFrames
import MAT
import StatsPlots
import NCDatasets

using DataFrames, StatsPlots, NCDatasets

include("base_GLODAP.jl")

# pathname = joinpath(@__DIR__, "woa18_decav_t00_01.nc")  # TEMP
# pathname = joinpath(@__DIR__, "woa18_decav_t13_01.nc")  # TEMP winter
pathname = joinpath(@__DIR__, "woa18_all_p00_01.nc")  # P
# pathname = joinpath(@__DIR__, "woa18_all_o00_01.nc")  # O2

ds = NCDataset(pathname) # get a list of global attributes

gr(size=(1000, 600))

# heatmap(ds["lon"][:], ds["lat"][:], ds["o_an"][:,:,2,1]')

contourf(ds["lon"][:], ds["lat"][:], ds["p_an"][:,:,1,1]', clims=(0,3),)

# savefig("WOA_TEMP.png")



