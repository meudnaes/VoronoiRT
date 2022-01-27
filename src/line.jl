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
    λ::Vector{Unitful.Quantity{T, Unitful.𝐋}}
    λidx::Vector{Int}
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
                            nλ_bb::Int, nλ_bf::Int,
                            gu::Int, gl::Int, f_value::T, atom_weight::Unitful.Mass{T},
                            Z::Int, temperature::Array{<: Unitful.Temperature})  where T <: AbstractFloat
        # Add conversion from cm^-1 to aJ, if type of χu is L^-1
        χu = Transparency.wavenumber_to_energy(χu)
        χl = Transparency.wavenumber_to_energy(χl)
        χ∞ = Transparency.wavenumber_to_energy(χ∞)
        @assert χ∞ > χu
        @assert χu > χl
        @assert gu > 0
        @assert gl > 0
        @assert f_value > 0
        @assert atom_weight > 0u"kg"
        @assert Z >= 1
        # Sample wavelengths for bound-bound and bound-free transitions
        λ0 = convert(Quantity{T, Unitful.𝐋}, ((h * c_0) / (χu - χl)) |> u"nm")
        qwing = 600.0
        qcore = 15.0
        λbb = sample_λ_line(nλ_bb, λ0, qwing, qcore)
        # from Ida
        ## Same χ for both
        λ1_min = transition_λ(χl, χ∞)*(1/2.0)^2 .+ 0.001u"nm"
        λ2_min = transition_λ(χl, χ∞)*(2/2.0)^2 .+ 0.001u"nm"
        ##
        λbf_l = sample_λ_boundfree(nλ_bf, λ1_min, χl, χ∞)
        λbf_u = sample_λ_boundfree(nλ_bf, λ2_min, χu, χ∞)
        λ = vcat(λbb, λbf_l, λbf_u)
        λi = [1, nλ_bb+1, nλ_bb+nλ_bf+1, nλ_bb+2*nλ_bf+1]
        # Einstein coefficients
        Aul = convert(Quantity{T, Unitful.𝐓^-1}, calc_Aji(λ0, gl / gu, f_value))
        Bul = calc_Bji(λ0, Aul)
        Blu = gu / gl * Bul
        # Doppler doppler_width
        ΔD = doppler_width.(λ0, atom_weight, temperature)

        new{T}(Aul, Bul, Blu, λ0, λ, λi, χl, χu, χ∞, gl, gu, atom_weight, Z, ΔD)
    end
end

"""
Computes the Voigt line profile
"""
function compute_voigt_profile(line::HydrogenicLine, atmos::Atmosphere,
                               damping_λ::Array{Float64, 4},
                               θ::Float64, ϕ::Float64)

    k = [cos(θ), cos(ϕ)*sin(θ), sin(ϕ)*sin(θ)]

    # calculate line of sight velocity
    v_los = line_of_sight_velocity(atmos, k)
    v = Array{Float64, 4}(undef, (length(line.λ), size(v_los)...))
    for l in eachindex(line.λ)
        v[l, :, :, :] = (line.λ[l] .- line.λ0 .+ line.λ0.*v_los./c_0)./line.ΔD .|> Unitful.NoUnits
    end

    # calculate line profile
    profile = Array{Float64, 4}(undef, (length(line.λ), size(v_los)...))u"m^-1"
    for l in eachindex(line.λ)
        profile[l, :, :, :] = voigt_profile.(damping_λ[l, :, :, :], v[l, :, :, :], line.ΔD)
    end

    return profile
end
function compute_voigt_profile(line::HydrogenicLine, sites::VoronoiSites,
                               damping_λ::Array{Float64, 2},
                               θ::Float64, ϕ::Float64)

    k = [cos(θ), cos(ϕ)*sin(θ), sin(ϕ)*sin(θ)]

    # calculate line of sight velocity
    v_los = line_of_sight_velocity(sites, k)
    v = Array{Float64, 4}(undef, (length(line.λ), sites.n))
    for l in eachindex(line.λ)
        v[l, :] = (line.λ[l] .- line.λ0 .+ line.λ0.*v_los./c_0)./line.ΔD .|> Unitful.NoUnits
    end

    # calculate line profile
    profile = Array{Float64, 2}(undef, (length(line.λ), sites.n))u"m^-1"
    for l in eachindex(line.λ)
        profile[l, :] = voigt_profile.(damping_λ[l, :], v[l, :], line.ΔD)
    end

    return profile
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

    nλ_bb = 5
    nλ_bf = 10

    gl = 2
    gu = 8

    f_value = 4.162E-01

    atom_weight = 1.673557692882144e-27u"kg"

    Z = 1

    return χu, χl, χ∞, nλ_bb, nλ_bf, gu, gl, f_value, atom_weight, Z
end
