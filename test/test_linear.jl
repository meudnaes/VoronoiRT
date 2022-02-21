using Plots
using Interpolations
include("../src/functions.jl")

pyplot()

function test_bilinear()
    # Same as example on wikipedia
    # https://en.wikipedia.org/wiki/Bilinear_interpolation

    # Unit square
    x_bounds = (0, 1)
    y_bounds = (0, 1)

    #   [f(x1, y1)    f(x1, y2)
    #    f(x2, y1)    f(x2, y2)]

    values = [0. 1.
              1. 0.5]

    x = collect(LinRange(0,1,100))
    y = collect(LinRange(0,1,100))

    F = Matrix{Float64}(undef, (length(x), length(y)))

    for i in eachindex(x)
        for j in eachindex(y)
            F[i, j] = bilinear(x[i], y[j], x_bounds, y_bounds, values)
        end
    end

    heatmap(x, y, F, c=:jet, dpi=500, title="Test Bilinear",
            xaxis="x", yaxis="y", aspect_ratio=:equal)

    savefig("bilinear_test")
end

function test_trilinear()
    # Unit cube
    x_bounds = [0.0, 1.0]
    y_bounds = [0.0, 1.0]
    z_bounds = [0.0, 1.0]

    #   [[f(x1, y1)    f(x1, y2)
    #    f(x2, y1)    f(x2, y2)]

    values = Array{Float64, 3}(undef, (2, 2, 2))
    values[:, 1, 1] .= 0.0
    values[:, 1, 2] .= 1.0
    values[:, 2, 1] .= 1.0
    values[:, 2, 2] .= 0.5

    x = collect(LinRange(0.1, 0.99, 100))
    y = collect(LinRange(0.1, 0.99, 100))
    z = collect(LinRange(0.1, 0.99, 100))

    F = Array{Float64, 3}(undef, (length(x), length(y), length(z)))

    for k in eachindex(z)
        for i in eachindex(x)
            for j in eachindex(y)
                F[k, i, j] = trilinear(z[k], x[i], y[j], z_bounds, x_bounds, y_bounds, values)
            end
        end
    end

    for k in [1, 50, 100]
        heatmap(x, y, F[k, :, :], c=:jet, dpi=500, title="Test Trilinear $k",
                xaxis="x", yaxis="y")

        savefig("trilinear_test_$k")
    end
end

function trilinear(z_mrk::Float64, x_mrk::Float64, y_mrk::Float64,
                   z::Vector{Float64}, x::Vector{Float64}, y::Vector{Float64},
                   vals::AbstractArray)
    # Returns the index of the first value in a greater than or equal to x
    # Subtract by 1 to get coordinate of lower corner
    idz = searchsortedfirst(z, z_mrk) - 1
    idx = searchsortedfirst(x, x_mrk) - 1
    idy = searchsortedfirst(y, y_mrk) - 1

    # bounding corner coordinates
    z0 = z[idz]; z1 = z[idz+1]
    x0 = x[idx]; x1 = x[idx+1]
    y0 = y[idy]; y1 = y[idy+1]

    # difference between coordinates and interpolation point
    x_d = (x_mrk - x0)/(x1 - x0)
    y_d = (y_mrk - y0)/(y1 - y0)
    z_d = (z_mrk - z0)/(z1 - z0)

    # values at each corner (z is first index in data array)
    c000 = vals[idz, idx, idy]
    c010 = vals[idz, idx, idy+1]
    c100 = vals[idz, idx+1, idy]
    c110 = vals[idz, idx+1, idy+1]
    c001 = vals[idz+1, idx, idy]
    c011 = vals[idz+1, idx, idy+1]
    c101 = vals[idz+1, idx+1, idy]
    c111 = vals[idz+1, idx+1, idy+1]

    # interpolate in x direction
    c00 = c000*(1 - x_d) + c100*x_d
    c01 = c001*(1 - x_d) + c101*x_d
    c10 = c010*(1 - x_d) + c110*x_d
    c11 = c011*(1 - x_d) + c111*x_d

    # interpolate in y direction
    c0 = c00*(1 - y_d) + c10*y_d
    c1 = c01*(1 - y_d) + c11*y_d

    # intepolate in z direction
    c = c0*(1 - z_d) + c1*z_d
    return c
end

test_bilinear()
test_trilinear()
