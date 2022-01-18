using Plots
using Transparency
include("functions.jl")

"""
    LTE_populations(atom::Atom,
                    temperature::Array{<:Unitful.Temperature, 3},
                    electron_density::Array{<:NumberDensity, 3})
Given the atom density, calculate the atom populations according to LTE.
Tiago
"""
function LTE_populations(atom::Atom,
                         temperature::Array{<:Unitful.Temperature, 3},
                         electron_density::Array{<:NumberDensity, 3})
    χ = atom.χ
    g = atom.g
    atom_density = atom.density
    nz,nx,ny = size(atom_density)

    n_levels = length(χ)
    n_relative = ones(Float64, nz,nx,ny, n_levels)

    saha_const = (k_B / h) * (2π * m_e) / h
    saha_factor = 2 * ((saha_const * temperature).^(3/2) ./ electron_density) .|> u"m/m"

    for i=2:n_levels
        ΔE = χ[i] - χ[1]
        n_relative[:,:,:,i] = g[i] / g[1] * exp.(-ΔE ./ (k_B * temperature))
    end

    # Last level is ionised stage (H II)
    n_relative[:,:,:,n_levels] .*= saha_factor
    n_relative[:,:,:,1] = 1 ./ sum(n_relative, dims=4)[:,:,:,1]
    n_relative[:,:,:,2:end] .*= n_relative[:,:,:,1]

    return n_relative .* atom_density
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

line = AtomicLine(test_atom()...)

DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=4)...)

LTE_pops = LTE_populations(atom::Atom,
                           atmos.temperature,
                           atmos.electron_density)
