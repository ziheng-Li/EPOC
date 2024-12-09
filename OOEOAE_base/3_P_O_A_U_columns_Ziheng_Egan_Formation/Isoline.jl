module Isoline

function test_find_isolines(f = (x, y) -> sqrt((x-0.25)^2+(y-0.5)^2)-0.25)

    nx, ny = 10, 20

    coord_values = fill((NaN, NaN), ny, nx)
    for i in 1:ny 
        for j in 1:nx
            coord_values[i, j] = (j/nx, i/ny)
        end
    end

    return find_isolines(f, coord_values)
end

"""
    find_isolines(f, coord_values::Matrix{T}; length_tol=1e-6) -> isolines::Vector{T}[]

Find lines where `f(x::T...) == 0.0` and that cross a boundary in a 2D surface defined by `coord_values`

Searches for boundary crossings and then uses piecewise linear continuation
on the grid defined by `coord_values`, see https://en.wikipedia.org/wiki/Piecewise_linear_continuation, followed
by a 1D solve to find crossing points to accuracy length_tol.

Each element of `coord_values` should be a Tuple of length equal to the number of arguments expected by `f`, eg:
- a Tuple length 2 if `coord_values` represents a 2D surface, where 
  `f` has arguments `f(x1, x2)`.
- a Tuple length 3 if `coord_values` represents a 2D surface embedded in a 3D space, where 
  `f` has arguments `f(x1, x2, x3)`

`coord_values` may contain NaN to indicate a point not present in the surface
"""
function find_isolines(
    f, coord_values::Matrix{T};
    root_solve=ITP_solve, 
    length_tol=1e-6, 
    verbose=false, 
    maxiters=20
) where T

    isolines = Vector{T}[]
    boundary_simplexes = Simplex2D[]
    # find all boundary crossings
    for (ijbnd, ijinc, crsize, ijswap) in [
            (1,                     1,   size(coord_values, 2), identity), # boundary i=1 (top row of matrix)
            (size(coord_values, 1), -1,  size(coord_values, 2), identity), # boundary i = size(coord_values, 1) (bottom row of matrix)
            (1,                     1,   size(coord_values, 1), reverse), # boundary j=1 (leftmost column of matrix)
            (size(coord_values, 2), -1,  size(coord_values, 1), reverse), # boundary j = size(coord_values, 2) (rightmost column of matrix)
        ]

        append!(boundary_simplexes, _search_const_ij(f, coord_values,ijbnd, ijinc, crsize, ijswap))
    end

    # keep track of end of each isoline so we don't repeat it
    line_end_simplexes = Simplex2D[]

    for start_simplex in boundary_simplexes
        # check if this is the end of an isoline we've already done
        if !isnothing(findfirst(x -> simplexes_out_in_same(x, start_simplex), line_end_simplexes))
            continue
        end

         # println("start_simplex: ", start_simplex)

        current_simplex = start_simplex
        end_simplex = ZeroSimplex2D
       
        isoline_coords = T[]
        push!(isoline_coords, coord_values_f0_edge(f, coord_values, current_simplex, current_simplex.in_edge; root_solve, length_tol, verbose, maxiters))

        while current_simplex != ZeroSimplex2D
            push!(isoline_coords, coord_values_f0_edge(f, coord_values, current_simplex, current_simplex.out_edge; root_solve, length_tol, verbose, maxiters))
            end_simplex = current_simplex            
            current_simplex = _find_next_Simplex2D(f, coord_values, current_simplex)
        end

        push!(line_end_simplexes, end_simplex)
        push!(isolines, isoline_coords)
    end
   
    return isolines
end


"""
    interpolate_isoline(l::Vector{T}, n; n_start=1, n_end=n) -> li::Vector{T}
    interpolate_isoline(l::Vector{T}, n, fill_from_line; max_distance=Inf) -> li::Vector{T}

Interpolate so output line `li` has total of n points, optionally padding with NaN or filling with values supplied by `fill_from_line`

If `fill_from_line` is not supplied:
- 1:n_start-1 NaN
- n_start:n_end points interpolated from `l`, equally spaced according to Euclidian distance
- n_end+1:n NaN

If `fill_from_line` is supplied and is not `nothing`, then `n_start` and `n_end` are calculated based on the closest points in `fill_from_line`,
and 1:n_start and n_end:n are filled with values from `fill_from_line`
"""
function interpolate_isoline(l::Vector{T}, n; n_start=1, n_end=n) where T

    lambda = [0.0]
  
    n_l = length(l)
    for i in 2:n_l
        d = _dist_tuple(l[i], l[i-1])
        new_lambda = last(lambda) + d
        push!(lambda, new_lambda)
    end

    high_i_l = 2 # index of current upper point for interpolation
    li = T[]

    T_nan = ntuple(x->NaN, fieldcount(T))
    # pad with n_start - 1 NaN
    append!(li, fill(T_nan, n_start - 1))   
    push!(li, first(l))
    for i in (n_start+1):n_end
        target_lambda = (i-n_start)/(n_end-n_start)*last(lambda)

        while target_lambda > lambda[high_i_l] && high_i_l < n_l
            high_i_l += 1
        end

        low_i_l = high_i_l - 1
        d_n_l = lambda[high_i_l] - lambda[low_i_l]
        wt_low = (lambda[high_i_l] - target_lambda)/d_n_l
        wt_high = (target_lambda - lambda[low_i_l])/d_n_l
        push!(li, _add_tuple(wt_low .* l[low_i_l], wt_high .* l[high_i_l]))
    end
  
    # pad with n - n_end NaN
    append!(li, fill(T_nan, n - n_end))   

    return li
end

function interpolate_isoline(l::Vector{T}, n, fill_from_line; max_distance=Inf) where T
    if isnothing(fill_from_line)
        li = interpolate_isoline(l, n)
    else
        n_start = Isoline._find_nearest(fill_from_line, first(l); max_distance)
        n_end = Isoline._find_nearest(fill_from_line, last(l); max_distance)

        li = interpolate_isoline(l, n; n_start, n_end)

        if n_start > 1
            li[1:n_start] .= fill_from_line[1:n_start]
        end

        if n_end < n
            li[n_end:n] .= fill_from_line[n_end:n]
        end
    end

    return li
end

"""
    Find points where `f(x::T...) == 0.0` along line defined by the Vector of coordinates `l`
"""
function find_zeros_line(
    f, l::Vector{T}; 
    length_tol=1e-6, maxiters=20, verbose=false, root_solve=ITP_solve,
) where T
    f_zeros = T[]

    last_f = f(l[1]...)
    if iszero(last_f)
        push!(f_zeros, l[1])
    end

    for i in 2:length(l)
        next_f = f(l[i]...)
        if iszero(next_f)
            push!(f_zeros[l[i]])
        elseif last_f * next_f < 0
            c_zero = interp_f0_coords(
                f, (l[i-1], l[i]), (last_f, next_f);
                root_solve, length_tol, maxiters, verbose,
            )
            push!(f_zeros, c_zero)
        end
        last_f = next_f
    end

    return f_zeros
end


# find index of point in l that is nearest to p and closer than max_distance
# return nothing if no point found
function _find_nearest(l::Vector{T}, p::T; max_distance=Inf) where T
    current_nearest_i = nothing
    current_nearest_d = Inf
    for (i, lp) in enumerate(l)
        d = _dist_tuple(lp, p) 
        if d <= max_distance && d < current_nearest_d
            current_nearest_d = d
            current_nearest_i  = i
        end
    end

    return current_nearest_i
end

############################################################################################################
# Piecewise linear continuation in 2D to find isolines
# https://en.wikipedia.org/wiki/Piecewise_linear_continuation
###########################################################################################################

"""
    Simplex2D(coord_indices, f, in_edge, out_edge)
    Simplex2D((in_ci1, in_ci2), (in_f1, in_f2), other_point_ci::CartesianIndex{2}, other_point_f) -> s::Simplex2D

2D simplex (a triangle) with function values at corners and two edges labelled as in and out

`Simplex2D((in_ci1, in_ci2), (in_f1, in_f2), other_point_ci::CartesianIndex{2}, other_point_f)` constructs a Simplex2D
from coordinate indices and function values of in edge, and coord indices and function value of a third point.
"""
struct Simplex2D
    "coordinate indices for the three corners"
    coord_indices::NTuple{3, CartesianIndex{2}}
    "function values at corners"
    f_values::NTuple{3, Float64}
    "location in coord_indices list of in edge"
    in_edge::NTuple{2, Int64}
    "location in coord_indices list of out edge"
    out_edge::NTuple{2, Int64}
end

function Simplex2D((in_ci1, in_ci2), (in_f1, in_f2), other_point_ci::CartesianIndex{2}, other_point_f)
    # println("Simplex2D: (in_ci1, in_ci2) ", (in_ci1, in_ci2), "  (in_f1, in_f2) ", (in_f1, in_f2), " other_point_ci ", other_point_ci)
     
    # look at function value to get output edge
    if _safe_sign(in_f1) != _safe_sign(other_point_f)
        out_edge = (3, 1)
    elseif _safe_sign(in_f2) != _safe_sign(other_point_f)
        out_edge = (2, 3)
    end

    return Simplex2D((in_ci1, in_ci2, other_point_ci), (in_f1, in_f2, other_point_f), (1, 2), out_edge)
end

const ZeroSimplex2D = Simplex2D(ntuple(x->CartesianIndex(0, 0), 3), (NaN, NaN, NaN), (0, 0), (0, 0))

coord_indices_in_edge(s::Simplex2D) = (s.coord_indices[s.in_edge[1]], s.coord_indices[s.in_edge[2]])
coord_indices_out_edge(s::Simplex2D) = (s.coord_indices[s.out_edge[1]], s.coord_indices[s.out_edge[2]])

# compare out edge of s1 to in edge of s2
function simplexes_out_in_same(s1::Simplex2D, s2::Simplex2D)
    s1_out_ci = coord_indices_out_edge(s1)
    s2_in_ci = coord_indices_in_edge(s2)
    return (
        s1_out_ci[1] in s2_in_ci  &&
        s1_out_ci[2] in s2_in_ci 
    )
end

# find the other index in [1, 2, 3] given any two indexes
function point_not_in_edge(edge::NTuple{2, Int64})
    other_point_map = [0 3 2; 3 0 1; 2 1 0]
    return other_point_map[edge[1], edge[2]]
end

coord_ok(coord_values) = !any(isnan.(coord_values))
f_at_coord(f, coord_values) = coord_ok(coord_values) ? f(coord_values...) : NaN

# find coordinate values where f = 0 on specified edge
function coord_values_f0_edge(f, coord_values, s::Simplex2D, edge::Tuple{Int, Int}; kwargs...) 
    return interp_f0_coords(
        f, 
        (coord_values[s.coord_indices[edge[1]]], coord_values[s.coord_indices[edge[2]]]),
        (s.f_values[edge[1]], s.f_values[edge[2]]);
        kwargs...
    )
end

# Use 1D root solver root_solve to find the coordinate values where f(c...) = 0.0, along the 
# line between (c1, c2) with f values at end points (f1(c1...), f2(c2...))
function interp_f0_coords(
    f, (c1, c2), (f1, f2);
    root_solve, length_tol, maxiters, verbose
)
    p = (f, c1, c2)
    # function wrapper that takes parameter p (with f and coordinates of interval end points) and maps interval to 0 <= lambda <= 1
    function f_new(lambda, p)
        f, c1, c2 = p
        cnew = _add_tuple((1.0 - lambda).*c1, lambda.*c2)
        return f(cnew...)
    end
 
    # rescale tolerance from length to lambda (which has values 0, 1 at the end points)
    lambda_tol = length_tol / _dist_tuple(c1, c2)
    lambda, niter, converged = root_solve(f_new, (0.0, 1.0); p=p, y=(f1, f2), ϵ=lambda_tol, maxiters, verbose)
    cnew = _add_tuple((1.0 - lambda).*c1, lambda.*c2)

    verbose && @info "interp_f0_coords: cnew $cnew  (c1, c2) $((c1, c2))  niter $niter\n\n"

    !converged && @warn "interp_f0_coords: failed niter $niter"

    return cnew
end





# true for +ve x and zero x, false for -ve x
# _safe_sign(x) = sign(x) >= 0.0
_safe_sign(x) = (x >= 0.0)
# add two tuples element-by-element
_add_tuple(t1, t2) = Tuple(t1[i] + t2[i] for i in 1:length(t1))
# Euclidian distance between two tuples considered as points in a vector space
_dist_tuple(t1, t2) = sqrt(sum((t1[i] - t2[i])^2 for i in 1:length(t1)))


"""
    _search_const_ij(f, coord_values, ijbnd,  ijinc, crsize, ijswap) -> s::Vector{Simplex2D}

Search a row or column (eg a boundary) with index `ijbnd` for edges with a zero crossing of f 
and return a Vector of Simplex2D with these as the input edges

# Arguments
- `ijbnd`: index of row or column to search
- `ijinc`:  increment +/- 1 to add to find the adjacent interior row/column
- `crsize`: number of columns or rows
- `ijswap`: function to map (i, j) -> (i, j) (to search along a row) or (i, j) -> (j, i) (to search along a column)
"""
function _search_const_ij(f, coord_values, ijbnd::Integer,  ijinc::Integer, crsize::Integer, ijswap)
    start_simplexes = Simplex2D[]
    # search boundary i = ibnd, interior i = ibnd+iinc
    ci1 = CartesianIndex(ijswap((ijbnd, 1)))
    f1 = f_at_coord(f, coord_values[ci1])
    for ij=2:crsize
        ci2 = CartesianIndex(ijswap((ijbnd, ij)))
        f2 = f_at_coord(f, coord_values[ci2])
        if _safe_sign(f1) != _safe_sign(f2)
            interior_point_ci = CartesianIndex(ijswap((ijbnd+ijinc, ij)))
            interior_point_f = f_at_coord(f, coord_values[interior_point_ci])
            if coord_ok(coord_values[ci1]) && coord_ok(coord_values[ci2]) && coord_ok(coord_values[interior_point_ci])
                push!(start_simplexes, Simplex2D((ci1, ci2), (f1, f2), interior_point_ci, interior_point_f))
            end
        end
        ci1 = ci2
        f1 = f2
    end

    return start_simplexes
end

# find next Simplex by "reflecting" s::Simplex2D across out edge
function _find_next_Simplex2D(f, coord_values, s::Simplex2D)
 
    # println("_find_next_Simplex2D s: ", s)

    next_in_edge_ci = (s.coord_indices[s.out_edge[1]], s.coord_indices[s.out_edge[2]])
    # println("_find_next_Simplex2D next_in_edge_ci: ", next_in_edge_ci)
    next_in_edge_f = (s.f_values[s.out_edge[1]], s.f_values[s.out_edge[2]])
    prev_other_ci = s.coord_indices[point_not_in_edge(s.out_edge)]
    # println("_find_next_Simplex2D prev_other_ci: ", prev_other_ci)

    # three cases to try to find new index 
    if next_in_edge_ci[1][1] == next_in_edge_ci[2][1]
        new_ci = CartesianIndex(next_in_edge_ci[1][1] - (prev_other_ci[1] - next_in_edge_ci[1][1]), next_in_edge_ci[1][2])
    elseif next_in_edge_ci[1][2] == next_in_edge_ci[2][2]
        new_ci = CartesianIndex(next_in_edge_ci[1][1], next_in_edge_ci[1][2] - (prev_other_ci[2] - next_in_edge_ci[1][2]))
    else
        # diagonal
        inew = next_in_edge_ci[1][1] == prev_other_ci[1] ? next_in_edge_ci[2][1] : next_in_edge_ci[1][1]
        jnew = next_in_edge_ci[1][2] == prev_other_ci[2] ? next_in_edge_ci[2][2] : next_in_edge_ci[1][2]
        new_ci = CartesianIndex(inew, jnew)
    end

    if new_ci[1] in 1:size(coord_values, 1) && new_ci[2] in 1:size(coord_values, 2) && coord_ok(coord_values[new_ci])
        new_f = f(coord_values[new_ci]...)
        return Simplex2D(next_in_edge_ci, next_in_edge_f, new_ci, new_f)
    else
        return ZeroSimplex2D
    end
end

#########################################################################################
# 1D root finding
# These algorithms are available in Roots.jl and SimpleNonlinearSolve.jl, but:
# - the Roots.jl API doesn't seem to let you set termination condition of x tolerance
# - the SimpleNonlinearSolve.jl version (at least of FalsePosition) looks broken
# So the simplest route was to reimplement them
########################################################################################

"""
    ITP_solve(f, a, b; p=nothing, ϵ=eps(), y=(f(a, p), f(b, p)), κ₁=0.2/(b-a), κ₂=2, n₀=1, maxiters=64, verbose=false) -> (x, niter, converged)

Find the root x of the function `y(x) = f(x, p)` on the interval [a, b], given optional fixed parameters p,
tolerance on x of ϵ,  and hyperparameters κ₁, κ₂, n₀. The hyperparameters are 
subject to the following constraints:

0 < κ₁ < Inf             
1 ≤ κ₂ < 1 + (1 + √5)/2   
0 < n₀ < Inf             

A Julia implementation of the Interpolate, Truncate Project algorithm by TheLateKronos https://discourse.julialang.org/u/thelatekronos/summary
from https://discourse.julialang.org/t/julia-implementation-of-the-interpolate-truncate-project-itp-root-finding-algorithm/77739
as detailed in `https://en.wikipedia.org/wiki/ITP_method`.
Hyperparameter defaults are set on the basis of the reasoning provided
in `https://docs.rs/kurbo/0.8.1/kurbo/common/fn.solve_itp.html`
"""
function ITP_solve(
    f, a, b; 
    p=nothing,
    ϵ=eps(),
    y=(f(a, p), f(b, p)),
    κ₁=0.2/(b-a),
    κ₂=2,
    n₀=1,
    maxiters=64,
    verbose=false,
)
    0 < κ₁ < Inf             ||   error("κ₁ must be between 0 and ∞")
    1 ≤ κ₂ < 1 + (1 + √5)/2   ||   error("κ₂ must be between 1 and 1+ (1 + √5)/2 (1 + the golden ratio)")
    0 < n₀ < Inf             ||   error("n₀ must be between 0 and ∞")
    n_1div2 = ceil(Int, log2((b-a)/2ϵ))
    nₘₐₓ = n_1div2 + n₀

    b > a || error("interval must have b > a")
    y_a, y_b = y
    iszero(y_a) && return a, 0, true
    iszero(y_b) && return b, 0, true

    sign_yb = sign(y_b) 
    sign(y_a) == sign_yb  &&  error("sign(f(a)) = sign(f(b)). There is no guaranteed root in the given interval.")

    verbose && @info "ITP_solve niter 0   a $a   b $b  f(a) $y_a  f(b) $y_b"
  
    # regularize so y_a*sign_yb < 0, y_b*sign_yb > 0
    y_a *= sign_yb  # always +ve
    y_b *= sign_yb  # always -ve

    niter = 0

    while b-a > 2ϵ && niter <= maxiters
        # Calculating parameters:
        x_1div2 = (a+b)/2
        r = ϵ*2^(nₘₐₓ - niter) - (b-a)/2
        δ = κ₁*(b-a)^κ₂
        
        # Interpolation:
        x_f = /(y_b*a - y_a*b, y_b-y_a)
        
        # Truncation:
        σ = sign(x_1div2 - x_f)
        δ ≤ abs(x_1div2 - x_f) ? (x_t=x_f+σ*δ) : (x_t = x_1div2)
        
        # Projection:
        abs(x_t - x_1div2) ≤ r ? (x_ITP = x_t) : (x_ITP = x_1div2 - σ*r)
        
        # Updating Interval:
        y_ITP = f(x_ITP, p)*sign_yb
        niter += 1

        verbose && @info "ITP_solve niter $niter   a $a   b $b  c $x_ITP  f(c) $(y_ITP*sign_yb)"
        if y_ITP > 0
            b = x_ITP
            y_b = y_ITP
        elseif y_ITP < 0
            a = x_ITP
            y_a = y_ITP
        else
            a = b = x_ITP
        end
    end
    converged = (niter <= maxiters)
    return (a+b)/2, niter, converged
end

ITP_solve(f, (a, b); kwargs...) = ITP_solve(f, a, b; kwargs...)

"""
    false_positon_solve(f, a, b; p=nothing, y=(f(a, p), f(b, p)), ϵ, maxiters=100, verbose=false) -> (x, niter, converged)

Find the root x of the function `y(x) = f(x, p)` on the interval [a, b], given optional fixed parameters p,
tolerance on x of ϵ.
"""
function false_position_solve(f, x1, x2; p=nothing, y=(f(x1, p), f(x2, p)), ϵ, maxiters=100, verbose=false)

    f1, f2 = y

    (xnew, fnew, (x1, x2), (f1, f2), xlength, xnearest) =  false_position_update(f, p, (x1, x2), (f1, f2))
    niter = 1

    verbose && println("false_position_solve: niter $niter xnew $xnew fnew $fnew xlength $xlength xnearest $xnearest (x1, x2) $((x1, x2)) (f1, f2) $((f1, f2))")

    while niter <= maxiters && xnearest > ϵ && abs(fnew) > 0.0
        (xnew, fnew, (x1, x2), (f1, f2), xlength, xnearest) =  _false_position_update(f, p, (x1, x2), (f1, f2))
        niter += 1
        verbose && println("false_position_solve: niter $niter xnew $xnew fnew $fnew xlength $xlength xnearest $xnearest (x1, x2) $((x1, x2)) (f1, f2) $((f1, f2))")
    end

    converged = niter <= maxiters

    return xnew, niter, converged
end

false_position_solve(f, (x1, x2); kwargs...) = false_position_solve(f, x1, x2; kwargs...)

function _false_position_update(f, p, (x1, x2), (f1, f2))

    l = f2 - f1
    wt1 = f2/l
    wt2 = 1.0 - wt1
    xnew = wt1*x1 + wt2*x2

    xlength = abs(x2 - x1)
    xnearest = min(abs(x2-xnew), abs(x1-xnew))

    fnew = f(xnew, p)

    if _safe_sign(fnew) == _safe_sign(f2)
        # replace 2
        new_xbracket = (x1, xnew)
        new_fbracket = (f1, fnew)
    else
        # replace 1
        new_xbracket = (xnew, x2)
        new_fbracket = (fnew, f2)
    end

    return (xnew, fnew, new_xbracket, new_fbracket, xlength, xnearest)
end

end # module