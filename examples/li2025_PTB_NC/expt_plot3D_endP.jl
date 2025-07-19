# import MarchingCubes
# import GeometryBasics
import GLMakie

import PALEOboxes as PB
import PALEOmodel
import PALEOcopse

import XLSX
import DataFrames
import Interpolations as IL

######################################
# key variables controlling 3D plot appearance
# See https://blog.makie.org/blogposts/v0.20/ for Makie 0.20 updates

P_lims = (0.0, 5.0)
O_lims = (0.0, 3.5)
A_lims = (0.5, 3.5)

# P_lims = (0.0, 5.0)
# O_lims = (0.0, 2.0)
# A_lims = (0.5, 1.5)


# global variables we use to set plot appearance
glm_size=(500, 500)
# https://discourse.julialang.org/t/statically-rotate-3d-scene-in-makie-jl/74356
glm_azimuth = 4.2  # Makie default is 4.005530633326986
glm_elevation = 0.3 # Makie default is 0.39269908169872414

glm_linewidth = 1 
glm_eqb_markersize = 10  # marker size to use for eqb points

# Lighting, see https://docs.makie.org/dev/reference/scene/lighting/, defaults are:
# julia> axs.scene.lights
#   2-element Vector{Makie.AbstractLight}:
#   Makie.AmbientLight(Observable(RGB{Float32}(0.45f0,0.45f0,0.45f0)))
#   Makie.DirectionalLight(Observable(RGB{Float32}(0.5f0,0.5f0,0.5f0)), Observable(Float32[-0.45679495, -0.6293204, -0.6287243]), true)
# glm_lights = [
#     GLMakie.AmbientLight(GLMakie.RGBf(0.45,0.45,0.45)),
#     GLMakie.DirectionalLight(GLMakie.RGBf(0.5,0.5,0.5), GLMakie.Vec3f(-0.457, -0.629, -0.629), true),
# ]
# light from left, with light fixed relative to surface
glm_lights = [
    GLMakie.AmbientLight(GLMakie.RGBf(0.45,0.45,0.45)),
    GLMakie.DirectionalLight(GLMakie.RGBf(0.5,0.5,0.5), GLMakie.Vec3f(1.0, 0.5, -1.0), false),
]

################################################################

"""
    plot_3D_critical_manifold(dPdt_surf, A_grid, folds_dPdt_lines_POA; kwargs... ) -> fig, axs, plt

Create Makie fig, axs, plot and plot critical manifold (dP/dt=0 nullcline) + folds
"""
function plot_3D_critical_manifold(
    dPdt_surf, A_grid, folds_dPdt_lines_POA;
    A_lines=true,
    P_lims=P_lims,
    O_lims=O_lims,
    A_lims=A_lims,
    glm_size=glm_size,
    glm_azimuth=glm_azimuth,
    glm_elevation=glm_elevation,
    glm_linewidth=glm_linewidth,
    glm_lights=glm_lights,
)
    GLMakie.activate!(inline=false) # pop-out window

    # dPdt surface
    colormap = :thermal
    # # colormap = :plasma
    dPdt_surf_P, dPdt_surf_O, dPdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_surf)
    fig, axs, plt = GLMakie.surface(
        dPdt_surf_P, dPdt_surf_O, dPdt_surf_A;
        invert_normals=true,
        backlight=1.0, # multiplier for light calculation on backside of surface
        # z coord -> color
        color=dPdt_surf_A, 
        # colormap=(colormap, 0.9),  # (colormap, alpha) Makie 0.19
        colormap=(colormap, 0.7),  # (colormap, alpha) Makie 0.20
        colorrange=A_lims,
        # solid colour
        # color=fill((:blue, 0.4), size(dPdt_surf_P)...), # workaround for a Makie bug see https://discourse.julialang.org/t/colored-surface-in-makie/61009/4
        transparency=true,
        # shading=true, # deprecated in Makie 0.20
        shading=GLMakie.MultiLightShading, # Makie 0.20
        # diffuse=GLMakie.Vec3f(0.4), # Makie 0.19
        diffuse=GLMakie.Vec3f(1.5), # Makie 0.20
        # specular=GLMakie.Vec3f(0.4), # Makie 0.19
        specular=GLMakie.Vec3f(0.8), # Makie 0.20
        shininess=32.0,
        # Doesn't work?
        # lightposition = GLMakie.Vec3f(15, 0, 0), # Light comes from (0, 0, 15), i.e the sphere
        figure = (
            # resolution=(500, 500), # deprecated in Makie 0.20
            size=glm_size,
        ),
        axis = (
            type=GLMakie.Axis3,
            # perspectiveness=0.5,
            perspectiveness=0.25,
            # perspectiveness=0.0, 
            azimuth=glm_azimuth,
            elevation=glm_elevation,
            xlabel="P_norm", xlabelcolor=:blue,
            ylabel="O_norm", ylabelcolor=:red,
            zlabel="A_norm",
        )
    )

    GLMakie.xlims!(axs, P_lims...)
    GLMakie.ylims!(axs, O_lims...)
    GLMakie.zlims!(axs, A_lims...) 

    # NB: .= to update vector, not replace
    # also requires display(fig) to update
    axs.scene.lights .= glm_lights

    if A_lines
        # check surface is correct
        # draw the lines on the surface
        A_target_vals = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5] # to match grid ticks
        for A_t in A_target_vals
            k = findfirst(x -> x > A_t, A_grid)
            if !isnothing(k)
                Pbal_P, Pbal_O, A_n = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_surf[:, k])
                GLMakie.lines!(
                    axs, Pbal_P, Pbal_O, A_n.*ones(length(Pbal_P));
                    # z-coord -> colour
                    color=A_n.*ones(length(Pbal_P)), # z coord -> color
                    colormap=(colormap, 0.4),  # (colormap, alpha)
                    colorrange=A_lims,
                    # solid colour
                    # color=:black,
                    # color=(:blue, 0.4),
                    # linestyle=:dashdot,
                )
            end
        end
    end

    # folds in dP/dt
    for fold_dPdt in folds_dPdt_lines_POA
        GLMakie.lines!(
            axs, SolverFunctionsOOEOAE2.tuples_to_coords(fold_dPdt)...;
            color=:blue,
            linewidth=2*glm_linewidth,
            # linestyle=:dash,
        )
    end

    return fig, axs, plt
end


"""
    plot_3D_A_nullcline(axs, dAdt_surf)

Plot dA/dt = 0 surface
"""
function plot_3D_A_nullcline(
    axs, dAdt_surf
)

    dAdt_surf_P, dAdt_surf_O, dAdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dAdt_surf)
    GLMakie.surface!(
        axs, dAdt_surf_P, dAdt_surf_O, dAdt_surf_A;
        # color=(:black, 0.1), # (color, alpha)
        # color=fill((:black, 0.1), size(dAdt_surf_P)...),  # workaround for a Makie bug see https://discourse.julialang.org/t/colored-surface-in-makie/61009/4
        color=fill((:black, 0.2), size(dAdt_surf_P)...), # Makie 0.20 increase alpha
        transparency=true,
        # shading=true, # deprecated in Makie 0.20
        shading=GLMakie.MultiLightShading, # Makie 0.20
    )

    return nothing
end


"""
    plot_3D_O_nullcline(axs, dOdt_surf)

Plot dO/dt = 0 surface
"""
function plot_3D_O_nullcline(
    axs, dOdt_surf
)
    dOdt_surf_P, dOdt_surf_O, dOdt_surf_A = SolverFunctionsOOEOAE2.tuples_to_coords(dOdt_surf)
    GLMakie.surface!(
        axs, dOdt_surf_P, dOdt_surf_O, dOdt_surf_A;
        # (color, alpha)
        # color=fill((:green, 0.2), size(dOdt_surf_P)...),  # workaround for a Makie bug see https://discourse.julialang.org/t/colored-surface-in-makie/61009/4
        # color=fill((:red, 0.1), size(dOdt_surf_P)...),
        color=fill((:red, 0.2), size(dOdt_surf_P)...), # Makie 0.20 increase alpha
        transparency=true,
        # shading=true, # deprecated in Makie 0.20
        shading=GLMakie.MultiLightShading, # Makie 0.20
    )

end

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