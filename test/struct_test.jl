using Unitful

struct Atmosphere_test
    z::Vector{<:Unitful.Length}
    x::Vector{<:Unitful.Length}
    y::Vector{<:Unitful.Length}
    # temperature::Array{<:Unitful.Temperature, 3}
    # electron_density::Array{<:NumberDensity, 3}
    # hydrogen_populations::Array{<:NumberDensity, 3}
    nz::Integer
    nx::Integer
    ny::Integer
    Δx::Unitful.Length
    Δy::Unitful.Length
    function Atmosphere_test(z::Vector{<:Unitful.Length{T}},
                             x::Vector{<:Unitful.Length{T}},
                             y::Vector{<:Unitful.Length{T}}) where T <: AbstractFloat
                        # temperature::Array{<:Unitful.Temperature{T}, 3},
                        # electron_density::Array{<:NumberDensity{T}, 3},
                        # hydrogen_populations::Array{<:NumberDensity{T}, 3}) where T <: AbstractFloat
        #Funcy func
        nz = length(z)
        nx = length(x)
        ny = length(y)
        Δx = x[2] - x[1]
        Δy = y[2] - y[1]
        new{T}(z, x, y, nz, nx, ny, Δx, Δy)
               # temperature, electron_density, hydrogen_populations,
               # nz, nx, ny, Δx, Δy)
    end
end


z = collect(LinRange(0,1,11))u"m"
x = copy(z)
y = copy(z)

T = 100*rand(11, 3)u"K"

atmos = Atmosphere(z, x, y)
