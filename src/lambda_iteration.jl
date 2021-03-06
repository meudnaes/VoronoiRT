include("rates.jl")
include("radiation.jl")
include("broadening.jl")
include("populations.jl")
include("voronoi_utils.jl")
include("characteristics.jl")
include("irregular_ray_tracing.jl")

function J_λ_regular(S_λ::Array{<:UnitsIntensity_λ, 4},
                     α_cont::Array{<:PerLength, 3},
                     populations::Array{<:NumberDensity, 4},
                     atmos::Atmosphere,
                     line::HydrogenicLine,
                     quadrature::String)

    # Ω = (θ, φ), space angle
    weights, θ_array, ϕ_array, n_angles = read_quadrature(quadrature)

    J_λ = zero(S_λ)

    γ = γ_constant(line,
                   atmos.temperature,
                   (populations[:, :, :, 1] .+ populations[:, :, :, 2]),
                   atmos.electron_density)

    damping_λ = Array{Float64, 4}(undef, size(S_λ))
    Threads.@threads for l in eachindex(line.λ)
        damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    # total extinction
    # α_tot = Array{Float64, 4}(undef, size(damping_λ))u"m^-1"

    for i in 1:n_angles
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]

        profile = compute_voigt_profile(line, atmos, damping_λ, k)

        Threads.@threads for l in eachindex(line.λ)

            α_tot = αline_λ(line,
                            profile[l, :, :, :],
                            populations[:, :, :, 2],
                            populations[:, :, :, 1]) + α_cont

            if θ > 90
                I_0 =  B_λ.(line.λ[l], atmos.temperature[1,:,:])
                J_λ[l,:,:,:] .+= weights[i].*short_characteristics_up(k,
                                                                      S_λ[l,:,:,:],
                                                                      I_0,
                                                                      α_tot,
                                                                      atmos;
                                                                      n_sweeps=3)
            elseif θ < 90
                I_0 = zero(S_λ[l, 1, :, :])
                J_λ[l,:,:,:] .+= weights[i].*short_characteristics_down(k,
                                                                        S_λ[l,:,:,:],
                                                                        I_0,
                                                                        α_tot,
                                                                        atmos;
                                                                        n_sweeps=3)
            end
        end
    end

    return J_λ, damping_λ
end

function J_λ_voronoi(S_λ::Matrix{<:UnitsIntensity_λ},
                     α_cont::Vector{<:PerLength},
                     populations::Matrix{<:NumberDensity},
                     sites::VoronoiSites,
                     line::HydrogenicLine,
                     quadrature::String)

    # Ω = (θ, φ), space angle
    weights, θ_array, ϕ_array, n_points = read_quadrature(quadrature)

    J_λ = zero(S_λ)

    γ = γ_constant(line,
                   sites.temperature,
                   (populations[:, 1] .+ populations[:, 2]),
                   sites.electron_density)

    damping_λ = Matrix{Float64}(undef, size(S_λ))
    Threads.@threads for l in eachindex(line.λ)
       damping_λ[l, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    n_sweeps = 3

    for i in 1:n_points
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]

        profile = compute_voigt_profile(line, sites, damping_λ, k)

        Threads.@threads for l in eachindex(line.λ)

            α_tot = αline_λ(line,
                            profile[l, :],
                            populations[:, 2],
                            populations[:, 1]) + α_cont

            if θ_array[i] > 90
                bottom_layer = sites.layers_up[2] - 1
                bottom_layer_idx = sites.perm_up[1:bottom_layer]
                I_0 = B_λ.(line.λ[l], sites.temperature[bottom_layer_idx])
                J_λ[l,:] += weights[i]*Delaunay_upII(k, S_λ[l,:], I_0, α_tot, sites, n_sweeps)

            elseif θ_array[i] < 90
                top_layer = sites.layers_down[2] - 1
                I_0 = zeros(top_layer)u"kW*nm^-1*m^-2"
                J_λ[l,:] += weights[i]*Delaunay_downII(k, S_λ[l,:], I_0, α_tot, sites, n_sweeps)
            end

        end
    end
    return J_λ, damping_λ
end


function Λ_regular(ϵ::AbstractFloat,
                   maxiter::Integer,
                   atmos::Atmosphere,
                   line::HydrogenicLine,
                   quadrature::String,
                   DATA::String)

    # LTE populations
    LTE_pops = LTE_populations(line, atmos)

    # Initial populations
    populations = copy(LTE_pops)
    # zero_radiation_populations(line, atmos)

    # Find continuum extinction and absorption extinction (without Thomson and Rayleigh)
    α_cont = α_absorption.(line.λ0,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           LTE_pops[:,:,:,1].+LTE_pops[:,:,:,2],
                           LTE_pops[:,:,:,3]) .+
             α_scattering.(line.λ0,
                           atmos.electron_density,
                           LTE_pops[:,:,:,1])

    # destruction probability (Should I include line???)
    ελ = destruction(LTE_pops, atmos.electron_density, atmos.temperature, line)
    println("Minimum $(minimum(ελ)) destruction probability")

    # Start with the source function as the Planck function
    B_0 = Array{Float64, 4}(undef, (length(line.λ), size(α_cont)...))u"kW*m^-2*nm^-1"
    for l in eachindex(line.λ)
        B_0[l, :, :, :] = B_λ.(line.λ[l], atmos.temperature)
    end

    S_new = copy(B_0)

    S_old = zero(S_new)

    C = calculate_C(atmos, LTE_pops)

    i=0

    local J_new
    # check where ε < 5e-3, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter, DATA)
        #############################
        # Calculate radiation field #
        #############################
        S_old = copy(S_new)
        J_new, damping_λ = J_λ_regular(S_old, α_cont, populations,
                                       atmos, line, quadrature)

        for l in eachindex(line.λ)
            S_new[l,:,:,:] = (1 .- ελ).*J_new[l,:,:,:] .+ ελ.*B_0[l,:,:,:]
        end

        #############################
        #      Calculate rates      #
        #############################
        R = calculate_R(atmos, line, J_new, damping_λ, LTE_pops)

        #############################
        #    Update populations     #
        #############################
        populations = get_revised_populations(R, C, atmos.hydrogen_populations*1.0)

        ############################
        #Write current data to file#
        ############################
        write_to_file(populations[:, 2:end-1, 2:end-1, :], DATA)
        write_to_file(S_new[:, :, 2:end-1, 2:end-1], DATA)

        ############################
        #     Iteration is done    #
        ############################
        i+=1
    end

    if i == maxiter
        println("Did not converge inside scope")
        return J_new, S_new, α_cont, populations
    end

    println("Converged in $i iterations")

    return J_new, S_new, α_cont, populations
end

function Λ_voronoi(ϵ::AbstractFloat,
                   maxiter::Integer,
                   sites::VoronoiSites,
                   line::HydrogenicLine,
                   quadrature::String,
                   DATA::String)
    println("---Iterating---")

    # Start in LTE
    LTE_pops = LTE_populations(line, sites)

    # Initial populations
    populations = copy(LTE_pops)
    # zero_radiation_populations(line, sites)

    # Find continuum extinction and absorption extinction (only with Thomson and Rayleigh)
    α_cont = α_absorption.(line.λ0,
                           sites.temperature,
                           sites.electron_density*1.0,
                           LTE_pops[:,1].+LTE_pops[:,2],
                           LTE_pops[:,3]) .+
             α_scattering.(line.λ0,
                           sites.electron_density,
                           LTE_pops[:,1])

    # Start with the source function as the Planck function
    B_0 = Array{Float64, 2}(undef, (length(line.λ), sites.n))u"kW*m^-2*nm^-1"

    for l in eachindex(line.λ)
        B_0[l, :] = B_λ.(line.λ[l], sites.temperature)
    end

    S_new = copy(B_0)

    S_old = zero(S_new)

    # destruction probability (Should I include line???)
    ελ = destruction(LTE_pops, sites.electron_density, sites.temperature, line)
    println("Minimum $(minimum(ελ)) destruction probability")

    C = calculate_C(sites, LTE_pops)

    i=0

    local J_new
    # check where ε < 1e-2, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter, DATA)
        #############################
        # Calculate radiation field #
        #############################
        S_old = copy(S_new)
        J_new, damping_λ = J_λ_voronoi(S_old, α_cont, populations,
                                       sites, line, quadrature)

        for l in eachindex(line.λ)
            S_new[l,:] = (1 .- ελ).*J_new[l,:] .+ ελ.*B_0[l,:]
        end

        #############################
        #      Calculate rates      #
        #############################
        R = calculate_R(sites, line, J_new, damping_λ, LTE_pops)

        #############################
        #    Update populations     #
        #############################
        populations = get_revised_populations(R, C, sites.hydrogen_populations)

        ############################
        #Write current data to file#
        ############################
        write_to_file(populations, DATA)
        write_to_file(S_new, DATA)

        ############################
        #     Iteration is done    #
        ############################
        i+=1
    end

    if i == maxiter
        println("Did not converge inside scope")
        return J_new, S_new, α_cont, populations
    end

    println("Converged in $i iterations")

    return J_new, S_new, α_cont, populations
end

function criterion(S_new::Array{<:UnitsIntensity_λ, 4},
                   S_old::Array{<:UnitsIntensity_λ, 4},
                   ϵ::Float64,
                   i::Int,
                   maxiter::Int,
                   DATA::String)
    diff = 0
    nλ = size(S_new)[1]
    for l in 1:nλ
        l_diff = maximum(abs.(1 .- S_old[l,:,:,:]./S_new[l,:,:,:])) |> Unitful.NoUnits
        diff = max(diff, l_diff)
        if isnan(l_diff)
            println("NaN DIFF!, index $l")
        end
    end
    if i > 0
        println("   Rel. diff.: $diff")
    end
    println("Iteration $(i+1)...")

    # Write convergence history to file
    write_to_file(diff, i+1, DATA)

    diff > ϵ && i < maxiter
end

function criterion(S_new::Matrix{<:UnitsIntensity_λ},
                   S_old::Matrix{<:UnitsIntensity_λ},
                   ϵ::Float64,
                   i::Int,
                   maxiter::Int,
                   DATA::String)
    diff = 0
    nλ = size(S_new)[1]
    for l in 1:nλ
        l_diff = maximum(abs.(1 .- S_old[l, :]./S_new[l, :])) |> Unitful.NoUnits
        diff = max(diff, l_diff)
        if isnan(l_diff)
            println("NaN DIFF!, index $l")
        end
    end
    if i > 0
        println("   Rel. diff.: $diff")
    end
    println("Iteration $(i+1)...")

    # Write convergence history to file
    write_to_file(diff, i+1, DATA)

    diff > ϵ && i < maxiter
end
