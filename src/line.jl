using Test
using Plots
using Transparency
using LinearAlgebra

include("functions.jl")
include("voronoi_utils.jl")

"""
    HydrogenicLine{T <: AbstractFloat}

Structure containing all atomic information
"""
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
                            Z::Int, temperature::Array{<:Unitful.Temperature})  where T <: AbstractFloat
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
        nλ_bb = length(λbb)
        # from Ida
        ## Same χ for both
        λ1_min = transition_λ(χl, χ∞)*(1/2.0)^2 .+ 0.001u"nm"
        λ2_min = transition_λ(χl, χ∞)*(2/2.0)^2 .+ 0.001u"nm"
        ##
        λbf_l = sample_λ_boundfree(nλ_bf, λ1_min, χl, χ∞)
        λbf_u = sample_λ_boundfree(nλ_bf, λ2_min, χu, χ∞)
        λ = vcat(λbb, λbf_l, λbf_u)#, 500.0u"nm")
        # Account for 500 nm
        λi = [0, nλ_bb, nλ_bb+nλ_bf, nλ_bb+2*nλ_bf]#+1]
        # Einstein coefficients
        Aul = convert(Quantity{T, Unitful.𝐓^-1}, calc_Aji(λ0, gl / gu, f_value))
        Bul = calc_Bji(λ0, Aul)
        Blu = gu / gl * Bul
        # Doppler doppler_width
        ΔD = doppler_width.(λ0, atom_weight, temperature)
        @test all( Inf .> ustrip.(ΔD) .>= 0.0 )

        new{T}(Aul, Bul, Blu, λ0, λ, λi, χl, χu, χ∞, gl, gu, atom_weight, Z, ΔD)
    end
end

"""
    compute_voigt_profile(line::HydrogenicLine, atmos::Atmosphere,
                          damping_λ::Array{Float64, 4}, k::Vector{Float64})

Computes the Voigt line profile on a regular grid
"""
function compute_voigt_profile(line::HydrogenicLine, atmos::Atmosphere,
                               damping_λ::Array{Float64, 4}, k::Vector{Float64})

    # calculate line of sight velocity
    # Remember to use -k!, since k is moving towards the ray
    v_los = line_of_sight_velocity(atmos, -k)

    # calculate line profile.
    profile = Array{Float64, 4}(undef, (length(line.λ), size(v_los)...))u"m^-1"

    Threads.@threads for l in eachindex(line.λ)
        v = (line.λ[l] .- line.λ0 .+ line.λ0.*v_los./c_0)./line.ΔD .|> Unitful.NoUnits
        profile[l, :, :, :] = voigt_profile.(damping_λ[l, :, :, :], v, line.ΔD)
    end

    return profile
end

function compute_voigt_profile(line::HydrogenicLine, atmos::Atmosphere,
                               damping_λ::Array{Float64, 3}, k::Vector{Float64},
                               l::Int)

    # calculate line of sight velocity
    # Remember to use -k!, since k is moving towards the ray
    v_los = line_of_sight_velocity(atmos, -k)

    # calculate line profile.
    profile = Array{Float64, 3}(undef, size(v_los))u"m^-1"

    v = (line.λ[l] .- line.λ0 .+ line.λ0.*v_los./c_0)./line.ΔD .|> Unitful.NoUnits
    profile = voigt_profile.(damping_λ, v, line.ΔD)

    return profile
end

"""
    compute_voigt_profile(line::HydrogenicLine, sites::VoronoiSites,
                          damping_λ::Matrix{Float64}, k::Vector{Float64})

Computes the Voigt line profile on an irregular grid
"""
function compute_voigt_profile(line::HydrogenicLine, sites::VoronoiSites,
                               damping_λ::Matrix{Float64}, k::Vector{Float64})

    # calculate line of sight velocity
    # Remember to use -k!, since k is moving towards the ray
    v_los = line_of_sight_velocity(sites, -k)

    # calculate line profile
    profile = Array{Float64, 2}(undef, (length(line.λ), sites.n))u"m^-1"

    Threads.@threads for l in eachindex(line.λ)
        v = @. (line.λ[l] - line.λ0 + line.λ0*v_los/c_0)/line.ΔD .|> Unitful.NoUnits
        profile[l, :] = voigt_profile.(damping_λ[l, :], v, line.ΔD)
    end

    return profile
end

"""
    compute_doppler_profile(line::HydrogenicLine,
                                 atmos::Atmosphere,
                                 k::Vector{Float64})

Computes pure Doppler/thermal broadening. Used for testing.
"""
function compute_doppler_profile(line::HydrogenicLine,
                                 atmos::Atmosphere,
                                 k::Vector{Float64})

    # calculate line of sight velocity
    # Remember to use -k!, since k is moving towards the ray
    v_los = line_of_sight_velocity(atmos, -k)

    # calculate line profile.
    profile = Array{Float64, 4}(undef, (length(line.λ), size(v_los)...))u"m^-1"
    for l in eachindex(line.λ)
        #v = (line.λ[l] .- line.λ0 .+ line.λ0.*v_los./c_0)./line.ΔD .|> Unitful.NoUnits
        Δλ = line.λ[l] - line.λ0 .+ line.λ0.*v_los./c_0
        profile[l, :, :, :] = doppler_profile.(Δλ, line.ΔD)
    end

    return profile
end

function doppler_profile(Δλ::Unitful.Length, ΔλD::Unitful.Length)
    1/(sqrt(π)*ΔλD)*exp(-(Δλ/ΔλD)^2)
end

"""
    line_of_sight_velocity(atmos::Atmosphere, k::Vector)

Computes the line of sight velocity in all locations of the regular grid, given
the direction of the ray, k.
"""
function line_of_sight_velocity(atmos::Atmosphere, k::Vector{Float64})
    v_los = Array{Float64, 3}(undef, size(atmos.velocity_z))u"m/s"

    for jj in 1:length(atmos.y)
        for ii in 1:length(atmos.x)
            for kk in 1:length(atmos.z)
                velocity = [atmos.velocity_z[kk, ii, jj],
                            atmos.velocity_x[kk, ii, jj],
                            atmos.velocity_y[kk, ii, jj]]

                v_los[kk, ii, jj] = dot(velocity, k)
            end
        end
    end
    return v_los
end

"""
    line_of_sight_velocity(sites::VoronoiSites, k::Vector)

Computes the line of sight velocity in all locations of the irregular grid,
given the direction of the ray, k.
"""
function line_of_sight_velocity(sites::VoronoiSites, k::Vector{Float64})
    v_los = Vector{Float64}(undef, sites.n)u"m/s"

    for ii in 1:sites.n
        velocity = [sites.velocity_z[ii],
                    sites.velocity_x[ii],
                    sites.velocity_y[ii]]
        v_los[ii] = dot(velocity, k)
    end
    return v_los
end

"""
    αline_λ(line::HydrogenicLine,
            profile::Array{<:PerLength},
            n_j::Array{<:NumberDensity},
            n_i::Array{<:NumberDensity})

Compute line extinction given an `AtomicLine` struct, `profile` defined per wavelength,
and upper and lower population densities `n_j` and `n_i`.
"""
function αline_λ(line::HydrogenicLine,
                 profile::Array{<:PerLength},
                 n_j::Array{<:NumberDensity},
                 n_i::Array{<:NumberDensity})

    return @. (h*c_0/(4*π*line.λ0) * profile * (n_i * line.Bij - n_j * line.Bji)) |> u"m^-1"
end

"""
    test_atom(nλ_bb::Int, nλ_bf::Int)

Parameters for two-level atom
"""
function test_atom(nλ_bb::Int, nλ_bf::Int)
    χl = 0.0u"cm^-1"
    χu = 82258.211u"cm^-1"
    χ∞ = 109677.617u"cm^-1"

    gl = 2
    gu = 8

    f_value = 4.162e-1

    atom_weight = mass_H

    Z = 1

    return χu, χl, χ∞, nλ_bb, nλ_bf, gu, gl, f_value, atom_weight, Z
end

"""
    sample_λ_line(nλ::Int64, χl::Unitful.Energy, χu::Unitful.Energy,
                            qwing::Float64, qcore::Float64)


Get sampling wavelengths. Bound free wavelengths are
linearly sampled, while the bound-bound follow the
log-sampling from github.com/ITA-Solar/rh.
Taken from https://github.com/f0rmIdabel/SolarMCRT
"""
function sample_λ_line(nλ::Int64,
                       λ0::Unitful.Length,
                       qwing::Float64,
                       qcore::Float64)

    # Make sure odd # of bb wavelengths
    if nλ > 0 && nλ%2 == 0
        nλ += 1
    end

    # Either 1 or five or more wavelengths
    if 1 < nλ < 5
        nλ = 5
    end

    # Initialise wavelength array
    λ = Array{Float64,1}(undef, nλ)u"nm"

    # =================================================
    # Bound-bound transition
    # Follows github.com/ITA-Solar/rh/blob/master/getlambda.c
    # =================================================
    if nλ == 1
        λ[1] = λ0

    elseif nλ >= 5
        vmicro_char = 2.5u"km/s"

        n = nλ/2 # Questionable
        β = qwing/(2*qcore)
        y = β + sqrt(β*β + (β - 1.0)*n + 2.0 - 3.0*β)
        b = 2.0*log(y) / (n - 1)
        a = qwing / (n - 2.0 + y*y)

        center = (nλ÷2) + 1
        λ[center] = λ0
        q_to_λ = λ[center] * vmicro_char / c_0

        for w=1:(nλ÷2)
            Δλ = a*(w + (exp(b*w) - 1.0)) * q_to_λ
            λ[center-w] = λ[center] - Δλ
            λ[center+w] = λ[center] + Δλ
        end
    end

    return λ
end

"""
    sample_λ(nλ_bb::Int64, nλ_bf::Int64,
             χl::Unitful.Energy, χu::Unitful.Energy, χ∞::Unitful.Energy)

Get sampling wavelengths. Bound free wavelengths are
linearly sampled, while the bound-bound follow the
log-sampling from github.com/ITA-Solar/rh.
Taken from https://github.com/f0rmIdabel/SolarMCRT
"""
function sample_λ_boundfree(nλ::Int64,
                            λ_min::Unitful.Length,
                            χl::Unitful.Energy,
                            χ∞::Unitful.Energy)


    λ_max  = transition_λ(χl, χ∞)

    # Initialise wavelength array
    λ = Array{Float64,1}(undef, nλ)u"nm"

    # =================================================
    # Bound-free transitions
    # Linear spacing
    # =================================================
    if nλ == 1

        λ[1] = λ_max

    elseif nλ > 1
        Δλ = (λ_max - λ_min)/(nλ-1)
        λ[1] = λ_min

        for w=2:nλ
            λ[w] = λ[w-1] + Δλ
        end
    end

    return λ
end


"""
    transition_λ(χ1::Unitful.Energy, χ2::Unitful.Energy)

Get the corresponding wavelength for
the energy difference between two levels.
"""
function transition_λ(χ1::Unitful.Energy, χ2::Unitful.Energy)
    ((h * c_0) / (χ2-χ1)) |> u"nm"
end


"""
    destruction(LTE_pops::Array{<:NumberDensity},
                electron_density::Array{<:NumberDensity},
                temperature::Array{<:Unitful.Temperature},
                line::HydrogenicLine)

Calculate the photon destruction probability ε, per eq (3.98) in Rutten, 2003
"""
function destruction(LTE_pops::Array{<:NumberDensity},
                     electron_density::Array{<:NumberDensity},
                     temperature::Array{<:Unitful.Temperature},
                     line::HydrogenicLine)
    A21 = line.Aji
    B21 = line.Bji
    C21 = Cij(2, 1, electron_density, temperature, LTE_pops)
    B_λ0 = B_λ.(line.λ0, temperature)
    ε_λ0 = @. C21/(C21 + A21 + B21*B_λ0)
end

"""
    source_line(atmos::Atmosphere, line::HydrogenicLine, populations::Matrix{<:NumberDensity})

Compute the line the source function, from eq. 2.72 in Rutten's notes
"""
function source_line(atmos::Atmosphere, line::HydrogenicLine, populations::Array{<:NumberDensity, 4})

    gl = 2
    gu = 8

    nl = populations[:,:,:,1]
    nu = populations[:,:,:,2]

    ratio = @. gu*nl/(gl*nu)

    S_λ = @. 2*h*c_0^2/line.λ0^5 * 1/(ratio - 1)

end
