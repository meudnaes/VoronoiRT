using Test
using Plots
using Transparency
using LinearAlgebra

include("functions.jl")
include("voronoi_utils.jl")

# TODO
# Explicitly state all functions from Transparency ?
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

    # println(trapz(V_v, profile[:, 4, 4, 4].*line.ΔD[4, 4, 4]) |> Unitful.NoUnits) # 1.0002645422865621
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

function compute_doppler_profile(line::HydrogenicLine, atmos::Atmosphere,
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
Compute line extinction given an `AtomicLine` struct, `profile` defined per wavelength,
and upper and lower population densities `n_j` and `n_i`.
"""
function αline_λ(line::HydrogenicLine,
                 profile::Array{<:PerLength},
                 n_j::Array{<:NumberDensity},
                 n_i::Array{<:NumberDensity})

    return @. (h*c_0/(4*π*line.λ0) * profile * (n_i * line.Bij - n_j * line.Bji)) |> u"m^-1"
end

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
    sample_from_destruction(atmos::Atmosphere)

Sample Voronoi sites by using destruction probability as probability density.
"""
function sample_from_destruction(atmos::Atmosphere, n_sites::Int)

    nλ_bb = 0
    nλ_bf = 0
    line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., atmos.temperature)
    LTE_pops = LTE_populations(line, atmos)
    ελ = destruction(LTE_pops, atmos.electron_density, atmos.temperature, line)

    positions = rejection_sampling(n_sites, atmos, ustrip.(ελ))
end

"""
    sample_from_extinction(atmos::Atmosphere,
                                λ0::Unitful.Length,
                                n_sites::Int)

Sample Voronoi sites by using continuum extinction as probalility density.
"""
function sample_from_extinction(atmos::Atmosphere,
                                λ0::Unitful.Length,
                                n_sites::Int)

    populations = LTE_ionisation(atmos)

    # Find continuum extinction and absorption extinction (without Thomson and Rayleigh)
    α_cont = α_absorption.(λ0,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           populations[:,:,:,1].+populations[:,:,:,2],
                           populations[:,:,:,3]) .+
             α_scattering.(λ0,
                           atmos.electron_density,
                           populations[:,:,:,1])

    positions = rejection_sampling(n_sites, atmos, log10.(ustrip.(α_cont)))
    return positions
end


function sample_from_total_extinction(atmos::Atmosphere,
                                      n_sites::Int)

    populations = LTE_ionisation(atmos)
    line = HydrogenicLine(test_atom(1, 1)..., atmos.temperature)

    γ = γ_constant(line,
                   atmos.temperature,
                   (populations[:, :, :, 1] .+ populations[:, :, :, 2]),
                   atmos.electron_density)

    # damping_λ = Array{Float64, 4}(undef, (1, size(atmos.temperature)...))
    damping_λ = damping.(γ, line.λ0, line.ΔD)

    # Straight up
    k = -[1.0, 0.0, 0.0]
    v_los = line_of_sight_velocity(atmos, -k)
    v = (line.λ0.*v_los./c_0)./line.ΔD .|> Unitful.NoUnits
    profile = voigt_profile.(damping_λ[:, :, :], v, line.ΔD)
    α_line = αline_λ(line,
                     profile[:, :, :],
                     populations[:, :, :, 2],
                     populations[:, :, :, 1])

    α_cont =  α_absorption.(line.λ0,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           populations[:,:,:,1].+populations[:,:,:,2],
                           populations[:,:,:,3]) .+
              α_scattering.(line.λ0,
                           atmos.electron_density,
                           populations[:,:,:,1])

    positions = rejection_sampling(n_sites, atmos, log10.(ustrip.(α_line.+α_cont)))
    return positions
end


"""
    sample_from_temp_gradient(atmos::Atmosphere,
                              λ0::Unitful.Length,
                              n_sites::Int)

Sample Voronoi sites by using temperature gradient in z-direction as probalility
density.
"""
function sample_from_temp_gradient(atmos::Atmosphere,
                                   n_sites::Int)

    temperature = ustrip.(atmos.temperature)
    z = ustrip.(atmos.z)

    temp_gradient = copy(temperature)

    temp_gradient[1,:,:] = @. (temperature[2,:,:] - temperature[1,:,:])/(z[2] - z[1])

    for k in 2:length(z)-1
        temp_gradient[k,:,:] = @. (temperature[k+1,:,:] - temperature[k,:,:])/(z[k+1] - z[k])
    end

    temp_gradient[end,:,:] = @. (temperature[end,:,:] - temperature[end-1,:,:])/(z[end] - z[end-1])

    positions = rejection_sampling(n_sites, atmos, abs.(temp_gradient))
    return positions
end

function LTE_ionisation(atmos::Atmosphere)

    χl = 0.0u"cm^-1"
    χu = 82258.211u"cm^-1"
    χ∞ = 109677.617u"cm^-1"

    χl = Transparency.wavenumber_to_energy(χl)
    χu = Transparency.wavenumber_to_energy(χu)
    χ∞ = Transparency.wavenumber_to_energy(χ∞)

    χ = [χl, χu, χ∞]
    # Ionised hydrogen -> g = 1
    g = [2, 8, 1]
    atom_density = atmos.hydrogen_populations
    nz, nx, ny = size(atom_density)

    n_levels = 3
    n_relative = ones(Float64, nz, nx, ny, n_levels)

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

function LTE_ionisation(sites::VoronoiSites)

    χl = 0.0u"cm^-1"
    χu = 82258.211u"cm^-1"
    χ∞ = 109677.617u"cm^-1"

    χl = Transparency.wavenumber_to_energy(χl)
    χu = Transparency.wavenumber_to_energy(χu)
    χ∞ = Transparency.wavenumber_to_energy(χ∞)

    χ = [χl, χu, χ∞]
    # Ionised hydrogen -> g = 1
    g = [2, 8, 1]
    atom_density = sites.hydrogen_populations
    n_sites = length(atom_density)

    n_levels = 3
    n_relative = ones(Float64, n_sites, n_levels)

    saha_const = (k_B / h) * (2π * m_e) / h
    saha_factor = 2 * ((saha_const * sites.temperature).^(3/2) ./ sites.electron_density) .|> u"m/m"

    for i=2:n_levels
        ΔE = χ[i] - χ[1]
        n_relative[:,i] = g[i] / g[1] * exp.(-ΔE ./ (k_B * sites.temperature))
    end

    # Last level is ionised stage (H II)
    n_relative[:,n_levels] .*= saha_factor
    n_relative[:,1] = 1 ./ sum(n_relative, dims=2)[:,1]
    n_relative[:,2:end] .*= n_relative[:,1]

    return n_relative .* atom_density
end

function destruction(LTE_pops::Array{<:NumberDensity},
                     electron_density::Array{<:NumberDensity},
                     temperature::Array{<:Unitful.Temperature},
                     line::HydrogenicLine)
    # destruction, eq (3.98) in Rutten, 2003
    A21 = line.Aji
    B21 = line.Bji
    C21 = Cij(2, 1, electron_density, temperature, LTE_pops)
    B_λ0 = B_λ.(line.λ0, temperature)
    ε_λ0 = @. C21/(C21 + A21 + B21*B_λ0)
end
