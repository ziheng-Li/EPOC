import MarchingCubes
import GeometryBasics

import PALEOboxes as PB
import PALEOmodel
# import PALEOreactions
import PALEOcopse

######################################
# key variables controlling 3D plot appearance
# See https://blog.makie.org/blogposts/v0.20/ for Makie 0.20 updates

P_lims = (0.0, 5.0)
# O_lims = (0.0, 2.5)
# O_lims = (0.0, 2.0)
O_lims = (0.0, 1.5)
# A_lims = (1.0, 5.0)
A_lims = (2.0, 5.0)

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
        colorrange=(2.0, 5.0),
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
        A_target_vals = [2.5, 3.0, 3.5, 4.0, 4.5] # to match grid ticks
        for A_t in A_target_vals
            k = findfirst(x -> x > A_t, A_grid)
            if !isnothing(k)
                Pbal_P, Pbal_O, A_n = SolverFunctionsOOEOAE2.tuples_to_coords(dPdt_surf[:, k])
                GLMakie.lines!(
                    axs, Pbal_P, Pbal_O, A_n.*ones(length(Pbal_P));
                    # z-coord -> colour
                    color=A_n.*ones(length(Pbal_P)), # z coord -> color
                    colormap=(colormap, 0.4),  # (colormap, alpha)
                    colorrange=(2.0, 5.0),
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