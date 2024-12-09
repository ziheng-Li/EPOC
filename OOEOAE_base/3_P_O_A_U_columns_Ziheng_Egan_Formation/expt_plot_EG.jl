import XLSX
import DataFrames
import Interpolations as IL

"""
    plot_data_C_U_O

    read data and plot, return two plots
"""
function plot_data_C_U_O(filename)

    df_C = DataFrames.DataFrame(XLSX.readtable(filename,"PTB-C"))
    df_U = DataFrames.DataFrame(XLSX.readtable(filename,"PTB-U"))
    df_TEMP = DataFrames.DataFrame(XLSX.readtable(filename,"PTB-TEMP"))

    plot_PTB_C = scatter(df_C.Age, df_C.d13Ccarb, color=:gray, alpha=0.01, label="δ¹³Ccarb")
    plot_PTB_U = scatter(df_U.Age, df_U.d238U, color=:gray, alpha=0.1, label = "δ²³⁸Ucarb")
    plot_PTB_TEMP = scatter(df_TEMP.Age, df_TEMP.d18O, color=:gray, alpha=0.1, label = "δ¹⁸Oapatite", ylims=(17, 23), yflip=true)

    return plot_PTB_C, plot_PTB_U, plot_PTB_TEMP, df_C, df_U, df_TEMP
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

Input the Dataframe, the program delete the matrix rows that duplicate the elements in the first column
"""
function unique_rows(A::DataFrames.DataFrame, col) 
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
    moving_window_average

    To keep the data density the same along the time,
    Set the range of the window, moving window to get the average
    Need to sort the dataframe by :Age, large to small
    
"""
function moving_window_average(
    data::DataFrames.DataFrame,        # data is a df, need both Age and isotopic data
    field::Symbol,
    window_range::Float64;  # unit Myr
    start_age = -1000e6,              # if not given, start form the start Age of data
) 
    
    start_age = max(start_age, minimum(data.Age))
    start_index = 1
    new_Age = []
    new_data = []
    V_in_window = []
    start_index = findfirst(x -> x > start_age, data.Age)

    for i in start_index:length(data.Age)
        if (data[i, :Age] > (start_age + window_range)) || i == length(data.Age)
            push!(new_Age, start_age + window_range/2)
            push!(new_data, sum(V_in_window)/length(V_in_window))
            V_in_window = []
            start_age += window_range
        end
        push!(V_in_window, data[i, field])
    end
    
    scatter(data.Age, data[:, field])
    display(plot!(new_Age, new_data))

    return new_Age, new_data
end

