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

    populations = LTE_populations(atmos)

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

"""
    sample_from_total_extinction(atmos::Atmosphere,
                                 n_sites::Int)

Sample Voronoi sites by using total extinction (for rays moving straight
upwards) as probalility density.
"""
function sample_from_total_extinction(atmos::Atmosphere,
                                      n_sites::Int)

    populations = LTE_populations(atmos)
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

"""
    sample_from_hydrogen_ionisation(atmos::Atmosphere,
                                n_sites::Int)

Sample Voronoi sites by using ionised hydrogen in LTE as probalility density.
"""
function sample_from_ionised_hydrogen(atmos::Atmosphere,
                                         n_sites::Int)

    populations = LTE_populations(atmos)

    positions = rejection_sampling(n_sites, atmos, log10.(ustrip.(populations[:,:,:,end])))
    return positions
end


"""
Integrated extinction over angles in quadrature, use this for sampling sites.
"""
function sample_from_avg_ext(DATA::String, n_sites::Int)

    weights, θ_array, ϕ_array, _ = read_quadrature("../quadratures/ul7n12.dat")

    atmos, S_λ, populations = read_quantities(DATA; periodic=false)

    line = HydrogenicLine(test_atom(50, 20)..., atmos.temperature)

    LTE_pops = LTE_populations(line, atmos)

    # Find continuum extinction
    α_cont = α_absorption.(line.λ0,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           LTE_pops[:,:,:,1].+LTE_pops[:,:,:,2],
                           LTE_pops[:,:,:,3]) .+
             α_scattering.(line.λ0,
                           atmos.electron_density,
                           LTE_pops[:,:,:,1])

    γ = γ_constant(line,
                   atmos.temperature,
                   (populations[:, :, :, 1].+populations[:, :, :, 2]),
                   atmos.electron_density)

    α_int = zeros(size(α_cont))u"m^-1"

    λ0idx = argmin(abs.(line.λ .- line.λ0))

    for i in eachindex(weights)
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        damping_λ = Array{Float64, 4}(undef, size(S_λ))
        Threads.@threads for l in eachindex(line.λ)
            damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
        end

        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        profile = compute_voigt_profile(line, atmos, damping_λ, k)

        α_tot = αline_λ(line,
                        profile[λ0idx, :, :, :],
                        populations[:, :, :, 2],
                        populations[:, :, :, 1]) .+ α_cont

        α_int .+= weights[i].*α_tot
    end
    positions = rejection_sampling(n_sites, atmos, log10.(ustrip.(α_int)))
    return positions
end

function sample_from_destruction(atmos::Atmosphere)

    nλ_bb = 0
    nλ_bf = 0
    line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., atmos.temperature)
    LTE_pops = LTE_populations(line, atmos)
    ελ = destruction(LTE_pops, atmos.electron_density, atmos.temperature, line)

    return ustrip.(ελ)
end

function sample_from_logNH_invT(atmos::Atmosphere, n_sites::Int)

    temperature = atmos.temperature
    N_H = atmos.hydrogen_populations

    quantity = log10.(ustrip.(N_H)) .* ustrip.(temperature).^(-2/5)

    positions = rejection_sampling(n_sites, atmos, quantity)
end

function sample_from_logNH_invT_rootv(atmos::Atmosphere, n_sites::Int)

    temperature = atmos.temperature
    N_H = atmos.hydrogen_populations
    v_x = atmos.velocity_x
    v_y = atmos.velocity_y
    v_z = atmos.velocity_z

    v_sqrd = ustrip.(v_x).^2 .+ ustrip.(v_y).^2 .+ ustrip.(v_z).^2

    quantity = log10.(ustrip.(N_H)) .* ustrip.(temperature).^(-2/5) .* v_sqrd.^(1/3)

    positions = rejection_sampling(n_sites, atmos, quantity)
end

function sample_from_invNH_invT(atmos::Atmosphere, n_sites::Int)
    temperature = atmos.temperature
    N_H = atmos.hydrogen_populations

    quantity = @. log10(ustrip(N_H))^(-2) * ustrip(temperature)^(-2/5)

    positions = rejection_sampling(n_sites, atmos, quantity)
end
