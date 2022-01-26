using Plots
using Transparency
include("functions.jl")
include("voronoi_utils.jl")

struct HydrogenicLine{T <: AbstractFloat}
    Aji::Unitful.Frequency{T}
    # Units of Bij/Bji defined for J_lambda
    Bji::Unitful.Quantity{T, Unitful.𝐋 * Unitful.𝐓^2 / Unitful.𝐌}
    Bij::Unitful.Quantity{T, Unitful.𝐋 * Unitful.𝐓^2 / Unitful.𝐌}
    λ0::Unitful.Length{T}
    λline::Vector{Unitful.Length}
    χi::Unitful.Energy{T}
    χj::Unitful.Energy{T}
    # Properties of atom, not line, but keeping here for now
    χ∞::Unitful.Energy{T}
    gi::Int
    gj::Int
    atom_weight::Unitful.Mass{T}
    Z::Int
    ΔD::Array{Unitful.Quantity{T, Unitful.𝐋}}
    function HydrogenicLine(χu::Quantity{T}, χl::Quantity{T}, χ∞::Quantity{T},
                            gu::Int, gl::Int, f_value::T, atom_weight::Unitful.Mass{T},
                            Z::Int, temperature)  where T <: AbstractFloat
        χu = Transparency.wavenumber_to_energy(χu)
        χl = Transparency.wavenumber_to_energy(χl)
        χ∞ = Transparency.wavenumber_to_energy(χ∞)
        # Add conversion from cm^-1 to aJ, if type of χu is L^-1
        @assert χ∞ > χu
        @assert χu > χl
        @assert gu > 0
        @assert gl > 0
        @assert f_value > 0
        @assert atom_weight > 0u"kg"
        @assert Z >= 1
        qwing = 600.0
        qcore = 15.0
        λ0 = convert(Quantity{T, Unitful.𝐋}, ((h * c_0) / (χu - χl)) |> u"nm")
        nλ = 11
        λline = sample_λ_line(nλ, λ0, qwing, qcore)
        Aul = convert(Quantity{T, Unitful.𝐓^-1}, calc_Aji(λ0, gl / gu, f_value))
        Bul = calc_Bji(λ0, Aul)
        Blu = gu / gl * Bul
        ΔD = doppler_width.(λ0, atom_weight, temperature)
        new{T}(Aul, Bul, Blu, λ0, λline, χl, χu, χ∞, gl, gu, atom_weight, Z, ΔD)
    end
end

"""
Compute line extinction given an `AtomicLine` struct, `profile` defined per wavelength,
and upper and lower population densities `n_u` and `n_l`.
"""
function αline_λ(line::HydrogenicLine,
                 profile::Array{<:PerLength},
                 n_l::Array{<:NumberDensity},
                 n_u::Array{<:NumberDensity})
    (h .* c_0 / (4 .* π .* line.λ0) .* profile .* (n_l .* line.Bij .- n_u .* line.Bji)) .|> u"m^-1"
end

function test_atom()
    χl = 0.0u"cm^-1"
    χu = 82258.211u"cm^-1"
    χ∞ = 109677.617u"cm^-1"

    gl = 2
    gu = 8

    f_value = 4.162E-01

    atom_weight = 1.673557692882144e-27u"kg"

    Z = 1

    return χu, χl, χ∞, gu, gl, f_value, atom_weight, Z
end
