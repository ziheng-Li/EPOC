import DataFrames

using DataFrames

"""
    read_GLODAP

    Input filename, return the Dataframe
"""
function read_GLODAP(filename)
    # Assume this file is in a top level folder containing subfolders with data compilations    
    pathname = joinpath(@__DIR__, "GLODAPv2.2020_"*filename*".mat") 
    @info "Read data from GLobal Ocean Data Analysis Project: $pathname"
    matdata = MAT.matread(pathname)

    column_names = collect(keys(matdata))
    @info "The keys (column names) of the data: $(column_names)"

    df = DataFrame()

    # Iterate through the keys in the loaded data
    for key in keys(matdata)
        # Get the value associated with the current key
        value = matdata[key]

        if key[1:2] != "G2"
            continue
        end # exclude "expocode" and "expocodeno"

        # Check if the value is an array
        if isa(value, Array)
            # If it is an array, add it as a column to the DataFrame
            df[!, key] = value[:]
        end
    end

    return df
end

"""
Data_filter

Do the data filter for the GLODAPv2, options to filter the bottomdepth, latitude, interpolated data,

Options to go through the data points and identify the sections

If min_section_oxygen != NaN, but a filter value (unit μM or μmol/kg), then only output the data from the section
with minimum water column oxygen conc less than this value
"""
function Data_filter(
    filenames,
    filter_bottomdepth::Float64=6000.0,
    filter_lat::String="all",
    filter_interpo::String="all",
    min_section_oxygen=NaN,
    )

    df = read_GLODAP(filenames[1])
    for i = 2:length(filenames)
        df_temp = read_GLODAP(filenames[i])
        df = vcat(df, df_temp)
    end

    if filter_interpo == "all"
    else 
        # filter ignore 'interpolated' and 'no value' rows 
        # df = df[df.G2phosphatef .== 2, :]
        df = df[df[!,"G2"*filter_interpo*"f"] .== 2, :]
    end

    if filter_lat == "all"
    elseif  filter_lat == "highlat"
        df = df[(df[!,"G2latitude"] .< -60) .|| (df[!,"G2latitude"] .> 60), :]
    elseif  filter_lat == "lowlat"
        df = df[(-60 .< df[!,"G2latitude"] .< 60), :]
    elseif  filter_lat == "tropical"
        df = df[(-15 .< df[!,"G2latitude"] .< 15), :]
    elseif  filter_lat == "subtropical"
        df = df[(15 .< df[!,"G2latitude"] .< 45) .|| (-45 .< df[!,"G2latitude"] .< -15), :]
    elseif  filter_lat == "subpolar"
        df = df[(45 .< df[!,"G2latitude"] .< 60) .|| (-60 .< df[!,"G2latitude"] .< -45), :]
    end

    df = df[df[!,"G2bottomdepth"] .< filter_bottomdepth,:]

    @info "The total number of the datapoint is: $(nrow(df))"

    if isnan(min_section_oxygen)
        @info "We do not identify sections"
    else
        (num_section, df) = identify_each_section(df)
        @info "The total number of sections is: $(num_section)"
        df = df[df[!,"section_minimum_oxygen"] .<= min_section_oxygen,:]
    end

    return df
end

"""
identify_each_section

Go through the dataset, generate two columns for the Dataframe

1. the index of the section

2. the minimum oxygen conc of this secion
"""
function identify_each_section(df)
    num_row = nrow(df)
    minimum_oxygen_each_section = df.G2oxygen[1] # unit μM or μmol/kg
    section_index = 1 # or the number of section
    minimum_oxygen = Float64[]
    num_points_each_section = 1
    section_num = Int64[]
    append!(section_num, section_index) # first data point belonging to the first section

    for i = 2:num_row
        if (df.G2depth[i] < df.G2depth[i-1]) && (df.G2depth[i] <= 100)
            # then df_shelf.G2depth[i] should be the start point of the next section
            for j in 1:num_points_each_section
                append!(minimum_oxygen, minimum_oxygen_each_section)
            end

            # reset everything
            section_index += 1
            minimum_oxygen_each_section = df.G2oxygen[i]
            num_points_each_section = 1
        else
            # we need data points within a section rise monotonically with depth, haven't been tested
            num_points_each_section += 1
            minimum_oxygen_each_section = min(minimum_oxygen_each_section, df.G2oxygen[i])

            if i == num_row
                for j in 1:num_points_each_section
                    append!(minimum_oxygen, minimum_oxygen_each_section)
                end
            end
        end
        append!(section_num, section_index)
    end
    df[!, "section_num"] =  section_num[:]
    df[!, "section_minimum_oxygen"] =  minimum_oxygen[:]

    return section_index, df
end

"""
plot_TEMP_O2_P_dC13

The rules of the input:

filenames: a vector of names of the ocean which should be a subvector of ["Arctic_Ocean" "Atlantic_Ocean" "Indian_Ocean" "Pacific_Ocean" "Merged_Master_File"]

filter_bottomdepth: only count the data from the sites that the maximum depth < filter_bottomdepth. e.g. filter_bottomdepth=200 m to represent `epicontinental` in GeoClim

filter_lat: only count the data from certain latitude, three options ["highlat", "lowlat", "all"]
    - "highlat" : lat > 60˚ || lat < -60˚
    - "lowlat"  : -60˚ < lat < 60˚
   
filter_interpo: use the columns have the postfix `f` to select [good_value, interpolated, no_value]
"""
function plot_TEMP_O2_P_dC13(
    filenames,
    filter_bottomdepth::Float64=6000.0,
    filter_lat::String="all",
    filter_interpo::String="all",
    )

    df = Data_filter(filenames, filter_bottomdepth, filter_lat, filter_interpo)

    gr(size=(1200, 800))

    p1 = histogram2d(df.G2temperature, df.G2depth, bins=(40, 40), yflip=true, title=string(filter_bottomdepth)*"(m)_"*filter_lat, xlabel="TEMP (degree C)", ylabel="Depth (m)", xlims=[0.0,30.0], ylims=[0.0,filter_bottomdepth], color=:plasma)

    p2 = histogram2d(df.G2oxygen, df.G2depth, bins=(40, 40), yflip=true, label=false, xlabel="Oxygen (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,400.0], ylims=[0.0,filter_bottomdepth], color=:plasma)

    p3 = histogram2d(df.G2phosphate, df.G2depth, bins=(40, 40), yflip=true, label=false, xlabel="Phosphate (mmol/m3)", ylabel="Depth (m)", xlims=[0.0,4.0], ylims=[0.0,filter_bottomdepth], color=:plasma)

    p4 = histogram2d(df.G2c13, df.G2depth, bins=(40, 40), yflip=true, label=false, xlabel="δ13C (‰)", ylabel="Depth (m)", xlims=[-5,5], ylims=[0.0,filter_bottomdepth], color=:plasma)

    l = @layout[a b; c d]

    plot(p1, p2, p3, p4, layout = l)

    # savefig(filename*"_only.png")

end