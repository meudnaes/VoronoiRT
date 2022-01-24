using Plots
using Transparency
include("functions.jl")

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
    function HydrogenicLine(χu::Quantity{T}, χl::Quantity{T}, χ∞::Quantity{T},
                        gu::Int, gl::Int, f_value::T, atom_weight::Unitful.Mass{T},
                        Z::Int)  where T <: AbstractFloat
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
        nλ = 41
        λline = sample_λ_line(nλ, λ0, qwing, qcore)
        Aul = convert(Quantity{T, Unitful.𝐓^-1}, calc_Aji(λ0, gl / gu, f_value))
        Bul = calc_Bji(λ0, Aul)
        Blu = gu / gl * Bul
        new{T}(Aul, Bul, Blu, λ0, λline, χl, χu, χ∞, gl, gu, atom_weight, Z)
    end
end

"""
    LTE_populations(atom::Atom,
                    temperature::Array{<:Unitful.Temperature, 3},
                    electron_density::Array{<:NumberDensity, 3})
Given the atom density, calculate the atom populations according to LTE.
Tiago
"""
function LTE_populations(line::HydrogenicLine,
                         atmos::Atmosphere)
    χ = [line.χi, line.χj, line.χ∞]
    # Ionised hydrogen -> g = 1
    g = [line.gi, line.gj, 1]
    atom_density = atmos.hydrogen_populations
    nz, nx, ny = size(atom_density)

    n_levels = 3
    n_relative = Array{Float64, 4}(undef, (nz, nx, ny, n_levels))

    saha_const = (k_B / h) * (2π * m_e) / h
    saha_factor = 2 * ((saha_const * atmos.temperature).^(3/2) ./ atmos.electron_density) .|> u"m/m"

    for i=2:n_levels
        ΔE = χ[i] - χ[1]
        n_relative[:,:,:,i] = g[i] / g[1] * exp.(-ΔE ./ (k_B * atmos.temperature))
    end

    # Last level is ionised stage (H II)
    n_relative[:,:,:,n_levels] .*= saha_factor
    n_relative[:,:,:,1] = 1 ./ sum(n_relative, dims=4)[:,:,:,1]
    n_relative[:,:,:,2:end] .*= n_relative[:,:,:,1]

    return n_relative .* atom_density
end

"""
Compute line extinction given an `AtomicLine` struct, `profile` defined per wavelength,
and upper and lower population densities `n_u` and `n_l`.
"""
function αline_λ(line::HydrogenicLine,
                 profile::Array{<:PerLength, 3},
                 n_u::Array{<:NumberDensity, 3},
                 n_l::Array{<:NumberDensity, 3})
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
