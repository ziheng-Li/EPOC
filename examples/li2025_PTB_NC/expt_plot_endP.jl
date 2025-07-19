import PALEOboxes as PB
import PALEOmodel
import PALEOcopse

import XLSX
import DataFrames
import Interpolations as IL

"""
    O2_U_to_fold_point

    obtain the position of the fold points from k_O2_U

    fold_point[1] = O2_U[end]
    fold_point[end] = O2_U[1] * CPsea_factor (e.g. CPsea_factor = 4)
"""
function O2_U_to_fold_point(O2_U, CPsea_factor)

    fold_point=[]

    pstart = O2_U[end]
    pend = O2_U[1]*CPsea_factor

    for i in 1:length(O2_U)
        push!(fold_point, pstart+(pend-pstart)*(i-1)/100)
    end

    return fold_point
end

"""
    plot_data_C_U

    read data and plot, return two plots and the isotope row_values

    The isotopic values are all sort by Age (rev=true) and unique_rows()
"""
function plot_data_C_U_O(filename)

    df_C = DataFrames.DataFrame(XLSX.readtable(filename,"PTB-C"));  sort!(df_C, :Age, rev=true)
    df_U = DataFrames.DataFrame(XLSX.readtable(filename,"PTB-U"));  sort!(df_U, :Age, rev=true)
    df_TEMP = DataFrames.DataFrame(XLSX.readtable(filename,"PTB-TEMP"));  sort!(df_TEMP, :Age, rev=true)

    plot_PTB_C = scatter(df_C.Age, df_C.d13Ccarb, color=:gray, alpha=0.01, label="δ¹³Ccarb")
    plot_PTB_U = scatter(df_U.Age, df_U.d238U, color=:gray, alpha=0.1, label = "δ²³⁸Ucarb")
    plot_PTB_TEMP = scatter(df_TEMP.Age, df_TEMP.TEMP, color=:gray, alpha=0.1, label = "δ¹⁸O_TEMP")

    return plot_PTB_C, plot_PTB_U, plot_PTB_TEMP, unique_rows(df_C, :Age), unique_rows(df_U, :Age), unique_rows(df_TEMP, :Age)
end

"""
    removemissing

    remove, missing and NaN elements from two columns, 
    once there is a missing or NaN elements, remove the row (both age and data)
"""
function removemissing(df, columns)
    index_v=[]
    length_df = size(df,1)
    for i in 1:length_df
        row_withmissing = false

        for col in columns
            # if !(isa(df[i, col], Float64)) & !(isa(df[i, col], Int64))
            #     @info "Input dataset has a unaccepted value, i=$(i), col=$(col), Type=$(typeof(df[i, col]))"
            # end
            if ismissing(df[i, col]) || isnan(df[i, col])
                row_withmissing = true
                break  # Exit the loop if a missing or NaN value is found in the current row
            end
        end

        if row_withmissing == false
            push!(index_v, i)
        end
    end

    return df[index_v,:]
end

"""
    unique_rows

Input the Dataframe, the program delete the matrix rows that duplicate the elements in the column col
"""
function unique_rows(A::DataFrames.DataFrame, col::Symbol) 
    # find the index of unique rows
    row_indices = Int[]
    row_values = []
    for i in 1:size(A, 1)
        v = A[i, col]
        if !(v in row_values)
            push!(row_indices, i)
            push!(row_values, v)
        end
    end

    return A[row_indices,:]
end

"""
    moving_average

    Two arguments: array and range (e.g. 21)
"""
function moving_average(data::Vector{T}, range::Int) where T
    pre_n = Int(floor(range/2))
    post_n = Int(ceil(range/2))

    moving_a = similar(data, T, length(data))

    for i in 1:length(data)
        # vector save the date in the range
        average_i = []

        for j in max(1, i-pre_n):min(length(data), i+post_n)
            push!(average_i, data[j])
        end

        moving_a[i] = sum(average_i)/length(average_i)
    end

    return moving_a
end

"""
    box_plot

    To keep the data density the same along the time,
    Set the range of the window, moving window to get the average
    Need to sort the dataframe by :Age, large to small

    plot mean, max, min of each box
    
"""
function box_plot(
    df_data::DataFrames.DataFrame,        # data is a df, need both Age and isotopic data
    field::Symbol,
    window_range::Float64;  # unit yr
    start_age = -1000e6,              # if not given, start form the start Age of data
    end_age = 0e6,
    ylabel = "",
    data_alpha = 0.1,
    yflip = false,
    ylims=(),
    marker=:+,
    color=:blue,
    plot_lines=false,
) 

    df_removemissing = deepcopy(removemissing(df_data, [:Age, field]))
    sort!(df_removemissing, :Age)
    df_uniquerow = unique_rows(df_removemissing, :Age)

    start_age = max(start_age, minimum(df_data.Age))
    end_age = min(end_age, maximum(df_data.Age))
    df_Age_filtered = df_uniquerow[start_age .< df_uniquerow.Age .< end_age, :]
    
    start_index = 1
    new_Age = []
    new_data_mean = []; new_data_max = []; new_data_min = []
    V_in_window = []
    start_index = findfirst(x -> x > start_age, df_Age_filtered.Age)

    for i in start_index:length(df_Age_filtered.Age)
        push!(V_in_window, df_Age_filtered[i, field])
        if (df_Age_filtered[i, :Age] > (start_age + window_range)) || i == length(df_Age_filtered.Age)
            push!(new_Age, start_age + window_range/2)
            push!(new_data_mean, sum(V_in_window)/length(V_in_window)) # 
            push!(new_data_min, minimum(V_in_window))
            push!(new_data_max, maximum(V_in_window))
            @info "Time = $(new_Age[end]), $(V_in_window))"
            V_in_window = []
            start_age += window_range
        end
    end

    new_Age = new_Age / -1e6
    data_age = df_Age_filtered.Age / -1e6

    gr(size=(500, 500)) 
    p1 = scatter(data_age, df_Age_filtered[:, field], color=color, alpha=data_alpha, ylims=ylims, ylabel=ylabel, xlabel="", marker=marker, yflip=yflip, label=false, legend=:topright, grid = false, xflip=true, xlims=[244,254], )
    if plot_lines
        plot!(p1, new_Age, new_data_mean, color=:red, linewidth=2, label=false, ) # label="box mean", 
        plot!(p1, new_Age, new_data_min, color=:black, alpha=1.0, linestyle = :dot, label=false) # label="box min", 
        plot!(p1, new_Age, new_data_max, color=:black, alpha=1.0, linestyle = :dot, label=false) # label="box max", 
    else
    end
    
    return p1, new_Age, new_data_mean, new_data_min, new_data_max
end


"""
    plot_dist_vs_

    Two arguments: 
        the Output_dist which save the output after the random search
        the symbol or field name we plot field vs dist
"""
function plot_dist_vs_(
    Output_dist::Vector{Any},
    field::Symbol;
    xlims=(-1,1),
    ylims=(0.0,1.0),
)
    field_values = [getfield(param, field) for param in Output_dist]
    dist_values = [param.dist for param in Output_dist]

    pp = scatter(field_values, dist_values, xlabel=field, ylabel="Distance", label=false)
    scatter!([getfield(Output_dist[1], field)],[Output_dist[1].dist], color=:red, 
        xlims = xlims, ylims=ylims, label=false)

    return pp
end