include("line.jl")
include("functions.jl")
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
                   (populations[:, :, :, 1].+populations[:, :, :, 2]),
                   atmos.electron_density)

    damping_λ = Array{Float64, 4}(undef, size(S_λ))
    for l in eachindex(line.λ)
        damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    for i in 1:n_angles
        θ = θ_array[i]
        ϕ = ϕ_array[i]

        profile = compute_voigt_profile(line, atmos, damping_λ, θ, ϕ)

        # total extinction
        α_tot = Array{Float64, 4}(undef, size(profile))u"m^-1"
        for l in eachindex(line.λ)
            α_tot[l,:,:,:] = αline_λ(line,
                                     profile[l, :, :, :],
                                     populations[:, :, :, 1],
                                     populations[:, :, :, 2])

            α_tot[l,:,:,:] += α_cont
        end

        for l in eachindex(line.λ)
            if θ_array[i] > 90
                I_0 =  B_λ.(line.λ[l], atmos.temperature[1,:,:])
                J_λ[l,:,:,:] .+= weights[i].*short_characteristics_up(θ_array[i],
                                                                      ϕ_array[i],
                                                                      S_λ[l,:,:,:],
                                                                      α_tot[l,:,:,:],
                                                                      atmos,
                                                                      degrees=true,
                                                                      I_0=I_0)
            elseif θ_array[i] < 90
                I_0 = zero(S_λ[l, 1, :, :])
                J_λ[l,:,:,:] .+= weights[i].*short_characteristics_down(θ_array[i],
                                                                        ϕ_array[i],
                                                                        S_λ[l,:,:,:],
                                                                        α_tot[l,:,:,:],
                                                                        atmos,
                                                                        degrees=true,
                                                                        I_0=I_0)
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
                   (populations[:, 1].+populations[:, 2]),
                   sites.electron_density)

    damping_λ = Matrix{Float64}(undef, size(S_λ))
    for l in eachindex(line.λ)
       damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    n_sweeps = 3

    for i in 1:n_points
        θ = θ_array[i]*π/180
        ϕ = ϕ_array[i]*π/180
        k = -[cos(θ), cos(ϕ)*sin(θ), sin(ϕ)*sin(θ)]

        profile = compute_voigt_profile(line, sites, damping_λ, θ, ϕ)

        # total extinction
        α_tot = Matrix{Float64}(undef, size(profile))u"m^-1"
        for l in eachindex(line.λ)
            α_tot[l,:] = αline_λ(line,
                                 profile[l, :],
                                 populations[:, 1],
                                 populations[:, 2])

            α_tot[l,:] += α_cont
        end

        for l in eachindex(line.λ)
            if θ_array[i] > 90
                bottom_layer = sites.layers_up[2] - 1
                bottom_layer_idx = sites.perm_up[1:bottom_layer]
                I_0 = B_λ.(500u"nm", sites.temperature[bottom_layer_idx])
                J_λ[l,:] += weights[i]*Delaunay_up(sites, I_0,
                                                   S_λ[l,:], α_tot[l,:], k, n_sweeps)
            elseif θ_array[i] < 90
                top_layer = sites.layers_down[2] - 1
                I_0 = zeros(top_layer)u"kW*nm^-1*m^-2"
                J_λ[l,:] += weights[i]*Delaunay_down(sites, I_0,
                                                     S_λ[l,:], α_tot[l,:], k, n_sweeps)
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

    # Start in LTE

    LTE_pops = LTE_populations(line, atmos)
    populations = copy(LTE_pops)

    # Find continuum extinction and absorption extinction (without Thomson and Rayleigh)
    α_cont = α_continuum.(line.λ0,
                          atmos.temperature*1.0,
                          atmos.electron_density*1.0,
                          populations[:, :, :, 1]*1.0,
                          populations[:, :, :, 3]*1.0)

    α_a = α_absorption.(line.λ0,
                        atmos.temperature*1.0,
                        atmos.electron_density*1.0,
                        populations[:, :, :, 1]*1.0,
                        populations[:, :, :, 3]*1.0)

    # destruction probability (Should I include line???)
    ελ = destruction(populations, atmos.electron_density, atmos.temperature, line)
    thick = ελ .> 1e-2

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
    while criterion(S_new, S_old, ϵ, i, maxiter, thick, DATA)
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
    populations = copy(LTE_pops)

    # Find continuum extinction and absorption extinction (only with Thomson and Rayleigh)
    α_cont = α_continuum.(line.λ0,
                          sites.temperature*1.0,
                          sites.electron_density*1.0,
                          populations[:, 1]*1.0,
                          populations[:, 3]*1.0)

    α_a = α_absorption.(line.λ0,
                        sites.temperature*1.0,
                        sites.electron_density*1.0,
                        populations[:, 1]*1.0,
                        populations[:, 3]*1.0)

    # Start with the source function as the Planck function
    B_0 = Array{Float64, 2}(undef, (length(line.λ), sites.n))u"kW*m^-2*nm^-1"

    for l in eachindex(line.λ)
        B_0[l, :] = B_λ.(line.λ[l], sites.temperature)
    end

    S_new = copy(B_0)

    S_old = zero(S_new)

    # destruction probability (Should I include line???)
    ελ = destruction(populations, sites.electron_density, sites.temperature, line)
    thick = ελ .> 1e-2

    C = calculate_C(sites, LTE_pops)

    i=0

    local J_new
    # check where ε < 1e-2, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter, thick, DATA)
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

        i+=1
    end

    if i == maxiter
        println("Did not converge inside scope")
        return J_new, S_new, α_cont, populations
    end

    println("Converged in $i iterations")

    return J_new, S_new, α_cont, populations
end

function destruction(LTE_pops::Array{<:NumberDensity},
                     electron_density::Array{<:NumberDensity},
                     temperature::Array{<:Unitful.Temperature},
                     line::HydrogenicLine)
    # destruction, eq (3.98) in Rutten, 2003
    A21 = line.Aji
    B21 = line.Bji
    C21 = Cij(2, 1, electron_density, temperature, LTE_pops)
    Bλ_0 = B_λ.(line.λ0, temperature)
    ελ_0 = C21./(C21 .+ A21 .+ B21.*Bλ_0)
end

function criterion(S_new::Array{<:UnitsIntensity_λ, 4},
                   S_old::Array{<:UnitsIntensity_λ, 4},
                   ϵ::Float64,
                   i::Int,
                   maxiter::Int,
                   indcs,
                   DATA::String)
    diff = 0
    nλ = size(S_new)[1]
    for l in 1:nλ
        l_diff = maximum(abs.((S_new[l, :, :, :][indcs] .- S_old[l, :, :, :][indcs])
                                ./S_new[l, :, :, :][indcs])) |> Unitful.NoUnits
        if l_diff > diff
            diff = l_diff
        end
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
                   indcs,
                   DATA::String)
    diff = 0
    nλ = size(S_new)[1]
    for l in 1:nλ
        l_diff = maximum(abs.((S_new[l, :][indcs] .- S_old[l, :][indcs])
                                ./S_new[l, :][indcs])) |> Unitful.NoUnits
        if l_diff > diff
            diff = l_diff
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
