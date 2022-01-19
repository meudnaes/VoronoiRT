using Plots
using Transparency
include("functions.jl")

struct HydrogenicLine{T <: AbstractFloat}
    Aji::Unitful.Frequency{T}
    # Units of Bij/Bji defined for J_lambda
    Bji::Unitful.Quantity{T, Unitful.𝐋 * Unitful.𝐓^2 / Unitful.𝐌}
    Bij::Unitful.Quantity{T, Unitful.𝐋 * Unitful.𝐓^2 / Unitful.𝐌}
    λ0::Unitful.Length{T}
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
        λ0 = convert(Quantity{T, Unitful.𝐋}, ((h * c_0) / (χu - χl)) |> u"nm")
        Aul = convert(Quantity{T, Unitful.𝐓^-1}, calc_Aji(λ0, gl / gu, f_value))
        Bul = calc_Bji(λ0, Aul)
        Blu = gu / gl * Bul
        new{T}(Aul, Bul, Blu, λ0, χl, χu, χ∞, gl, gu, atom_weight, Z)
    end
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

line = HydrogenicLine(test_atom()...)

DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=4)...)

LTE_pops = LTE_populations(line, atmos)
