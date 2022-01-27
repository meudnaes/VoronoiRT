include("atom.jl")
include("line.jl")
include("functions.jl")
include("broadening.jl")
include("populations.jl")
include("voronoi_utils.jl")
include("characteristics.jl")
include("irregular_ray_tracing.jl")

function J_λ_regular(S_λ::AbstractArray,
                     α_tot::AbstractArray,
                     atmos::Atmosphere,
                     quadrature::String)

    # Ω = (θ, φ), space angle
    weights, θ_array, ϕ_array, n_points = read_quadrature(quadrature)

    J = zero(S_λ)

    for i in 1:n_points
        if θ_array[i] > 90
            I_0 =  B_λ.(500u"nm", atmos.temperature[1,:,:])
            J += weights[i]*short_characteristics_up(θ_array[i], ϕ_array[i], S_λ,
                                                     α_tot, atmos, degrees=true, I_0=I_0)
        elseif θ_array[i] < 90
            I_0 = zero(S_λ[1, :, :])
            J += weights[i]*short_characteristics_down(θ_array[i], ϕ_array[i], S_λ,
                                                       α_tot, atmos, degrees=true, I_0=I_0)
        end
    end
    return J
end

function J_λ_regular(S_λ::Array{<:UnitsIntensity_λ, 4},
                     α_cont::Array{<:PerLength, 4},
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

        # line extinction
        α_line = Array{Float64, 4}(undef, size(profile))u"m^-1"
        for l in eachindex(line.λ)
            α_line[l, :, :, :] = αline_λ(line,
                                         profile[l, :, :, :],
                                         populations[:, :, :, 1],
                                         populations[:, :, :, 2])
        end

        # total exinction
        α_tot = α_line .+ α_cont

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

function J_λ_voronoi(S_λ::AbstractArray,
                     α_tot::AbstractArray,
                     sites::VoronoiSites,
                     quadrature::String)

    # Ω = (θ, φ), space angle
    weights, θ_array, ϕ_array, n_points = read_quadrature(quadrature)

    J = zero(S_λ)

    n_sweeps = 3
    Nran = 1

    for i in 1:n_points
        θ = θ_array[i]*π/180
        ϕ = ϕ_array[i]*π/180
        k = -[cos(θ), cos(ϕ)*sin(θ), sin(ϕ)*sin(θ)]
        if θ_array[i] > 90
            perm = sortperm(sites.layers_up)
            layers_sorted = sites.layers_up[perm]
            bottom_layer = searchsortedfirst(layers_sorted, 2)-1
            I_0 = B_λ.(500u"nm", sites.temperature[1:bottom_layer])
            J += weights[i]*Delaunay_up(sites, I_0,
                                           S_λ, α_tot, k, n_sweeps, Nran)
        elseif θ_array[i] < 90
            perm = sortperm(sites.layers_down)
            layers_sorted = sites.layers_down[perm]
            top_layer = searchsortedfirst(layers_sorted, 2)-1
            I_0 = zeros(top_layer)u"kW*nm^-1*m^-2"
            J += weights[i]*Delaunay_down(sites, I_0,
                                             S_λ, α_tot, k, n_sweeps, Nran)
        end
    end
    return J
end

function J_λ_voronoi(S_λ::Vector{<:UnitsIntensity_λ},
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

    damping_λ = Array{Float64, 4}(undef, size(S_λ))
    for l in eachindex(line.λ)
       damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    n_sweeps = 3

    for i in 1:n_points
        θ = θ_array[i]*π/180
        ϕ = ϕ_array[i]*π/180
        k = -[cos(θ), cos(ϕ)*sin(θ), sin(ϕ)*sin(θ)]

        profile = compute_voigt_profile(line, sites, damping_λ, θ, ϕ)

        # line extinction
        α_line = Matrix{PerLength}(undef, size(profile))
        for l in eachindex(line.λ)
            α_line[l, :] = αline_λ(line,
                                   profile[l, :],
                                   populations[:, 1],
                                   populations[:, 2])
        end

        α_tot = α_cont .+ α_line
        for l in eachindex(line.λ)
            if θ_array[i] > 90
                perm = sortperm(sites.layers_up)
                layers_sorted = sites.layers_up[perm]
                bottom_layer = searchsortedfirst(layers_sorted, 2)-1
                I_0 = B_λ.(500u"nm", sites.temperature[1:bottom_layer])
                J += weights[i]*Delaunay_up(sites, I_0,
                                            S_λ[l,:], α_tot[l,:], k, n_sweeps)
            elseif θ_array[i] < 90
                perm = sortperm(sites.layers_down)
                layers_sorted = sites.layers_down[perm]
                top_layer = searchsortedfirst(layers_sorted, 2)-1
                I_0 = zeros(top_layer)u"kW*nm^-1*m^-2"
                J += weights[i]*Delaunay_down(sites, I_0,
                                              S_λ[l,:], α_tot[l,:], k, n_sweeps)
            end
        end
    end
    return J
end

function Λ_regular(ϵ::AbstractFloat,
                   maxiter::Integer,
                   atmos::Atmosphere,
                   quadrature::String)
    # choose a wavelength
    λ = 500u"nm"  # nm

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_tot = α_continuum.(λ, atmos.temperature*1.0, atmos.electron_density*1.0,
                    atmos.hydrogen_populations*1.0, atmos.hydrogen_populations*1.0)

    α_a = α_absorption.(λ, atmos.temperature*1.0, atmos.electron_density*1.0,
                        atmos.hydrogen_populations*1.0, atmos.hydrogen_populations*1.0)

    # destruction
    ε_λ = α_a ./ α_tot

    thick = ε_λ .> 5e-3

    # Start with the source function as the Planck function
    B_0 = B_λ.(λ, atmos.temperature)
    S_new = B_0

    S_old = zero(S_new)

    i=0

    local J_new
    # check where ε < 1e-2, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter, thick)
        print("Iteration $(i+1)\r")
        S_old = copy(S_new)
        J_new = J_λ_regular(S_old, α_tot, atmos, quadrature)
        S_new = (1 .- ε_λ).*J_new .+ ε_λ.*B_0
        i+=1
    end

    if i == maxiter
        println("Did not converge inside scope")
        return J_new, S_new, α_tot
    end

    println("Converged in $i iterations")

    return J_new, S_new, α_tot
end

function Λ_regular(ϵ::AbstractFloat,
                   maxiter::Integer,
                   atmos::Atmosphere,
                   line::HydrogenicLine,
                   quadrature::String)

    # Start in LTE

    populations = LTE_populations(line, atmos)*1.0

    # Find continuum extinction and absorption extinction (only with Thomson and Rayleigh)
    α_cont = Array{Float64, 4}(undef, (length(line.λ), size(atmos.temperature)...))u"m^-1"
    α_a = copy(α_cont)
    for l in eachindex(line.λ)
        α_cont[l, :, :, :] = α_continuum.(line.λ[l],
                                          atmos.temperature*1.0,
                                          atmos.electron_density*1.0,
                                          populations[:, :, :, 1]*1.0,
                                          populations[:, :, :, 3]*1.0)

        α_a[l, :, :, :] = α_absorption.(line.λ[l],
                                        atmos.temperature*1.0,
                                        atmos.electron_density*1.0,
                                        populations[:, :, :, 1]*1.0,
                                        populations[:, :, :, 3]*1.0)
    end

    # destruction probability (Should I include line???)
    ελ = Array{Float64, 4}(undef, size(α_cont))
    for l in eachindex(line.λ)
        ελ[l, :, :, :] = destruction(populations, atmos.electron_density, atmos.temperature, line)
    end
    thick = ελ .> 1e-2

    # Start with the source function as the Planck function
    B_0 = Array{Float64, 4}(undef, size(α_cont))u"kW*m^-2*nm^-1"

    for l in eachindex(line.λ)
        B_0[l, :, :, :] = B_λ.(line.λ[l], atmos.temperature)
    end

    S_new = copy(B_0)

    S_old = zero(S_new)

    i=0

    local J_new
    # check where ε < 5e-3, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter, thick)
        #############################
        # Calculate radiation field #
        #############################
        S_old = copy(S_new)
        J_new, damping_λ = J_λ_regular(S_old, α_cont, populations,
                                       atmos, line, quadrature)
        S_new = (1 .- ελ).*J_new .+ ελ.*B_0

        #############################
        #      Calculate rates      #
        #############################
        R, C = calculate_transition_rates(atmos, line, J_new, damping_λ)

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
                   quadrature::String)
    println("---Iterating---")

    # Start in LTE
    populations = LTE_populations(line, sites)

    # Find continuum extinction and absorption extinction (only with Thomson and Rayleigh)
    α_cont = Array{Float64, 2}(undef, (length(line.λ), sites.n))u"m^-1"
    α_a = copy(α_cont)
    for l in eachindex(line.λ)
        α_cont[l, :] = α_continuum.(line.λ[l],
                                          sites.temperature*1.0,
                                          sites.electron_density*1.0,
                                          populations[:, 1]*1.0,
                                          populations[:, 3]*1.0)

        α_a[l, :] = α_absorption.(line.λ[l],
                                        sites.temperature*1.0,
                                        sites.electron_density*1.0,
                                        populations[:, 1]*1.0,
                                        populations[:, 3]*1.0)
    end

    # Start with the source function as the Planck function
    B_0 = Array{Float64, 2}(undef, size(α_cont))u"kW*m^-2*nm^-1"

    for l in eachindex(line.λ)
        B_0[l, :] = B_λ.(line.λ[l], sites.temperature)
    end

    S_new = copy(B_0)

    S_old = zero(S_new)

    # destruction probability (Should I include line???)
    ελ = Array{Float64, 2}(undef, size(α_cont))
    for l in eachindex(line.λ)
        ελ[l, :] = destruction(populations, sites.electron_density, sites.temperature, line)
    end
    thick = ελ .> 1e-2

    i=0

    local J_new
    # check where ε < 1e-2, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter, thick)
        #############################
        # Calculate radiation field #
        #############################
        S_old = copy(S_new)
        J_new, damping_λ = J_λ_voronoi(S_old, α_cont, populations,
                                       sites, line, quadrature)

        S_new = (1 .- ελ).*J_new .+ ελ.*B_0

        #############################
        #      Calculate rates      #
        #############################
        R, C = calculate_transition_rates(sites, line, J_new, damping_λ)

        #############################
        #    Update populations     #
        #############################
        populations = get_revised_populations(R, C, sites.hydrogen_populations)

        i+=1
    end

    if i == maxiter
        println("Did not converge inside scope")
        return J_new, S_new, α_tot, populations
    end

    println("Converged in $i iterations")

    return J_new, S_new, α_tot, populations
end

function Λ_voronoi(ϵ::AbstractFloat, maxiter::Integer, sites::VoronoiSites, quadrature::String)
    println("---Iterating---")
    # choose a wavelength
    λ = 500u"nm"  # nm

    # Only continuum
    η_ν = 0

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_tot = α_continuum.(λ, sites.temperature*1.0, sites.electron_density*1.0,
                    sites.hydrogen_populations*1.0, sites.hydrogen_populations*1.0)

    α_a = α_absorption.(λ, sites.temperature*1.0, sites.electron_density*1.0,
                        sites.hydrogen_populations*1.0, sites.hydrogen_populations*1.0)

    # destruction
    ε_λ = α_a ./ α_tot

    thick = ε_λ .> 5e-3

    # Start with the source function as the Planck function
    B_0 = B_λ.(λ, sites.temperature)
    S_new = B_0

    S_old = zero(S_new)

    i=0

    local J_new
    # check where ε < 1e-2, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter, thick)
        print("Iteration $(i+1)\r")
        S_old = copy(S_new)
        J_new = J_λ_voronoi(S_old, α_tot, sites, quadrature)
        S_new = (1 .- ε_λ).*J_new .+ ε_λ.*B_0
        i+=1
    end

    if i == maxiter
        println("Did not converge inside scope")
        return J_new, S_new, α_tot
    end

    println("Converged in $i iterations")

    return J_new, S_new, α_tot
end

function destruction(populations::Array{<:NumberDensity},
                     electron_density::Array{<:NumberDensity},
                     temperature::Array{<:Unitful.Temperature},
                     line::HydrogenicLine)
    # destruction, eq (3.98) in Rutten, 2003
    A21 = line.Aji
    B21 = line.Bji
    C21 = Cij(2, 1, electron_density, temperature, populations)
    Bλ_0 = B_λ.(line.λ0, temperature)
    ελ_0 = C21./(C21 .+ A21 .+ B21.*Bλ_0)
end

function criterion(S_new, S_old, ϵ, i, maxiter, indcs)
    diff = maximum(abs.((S_new[indcs] .- S_old[indcs])./S_new[indcs])) |> Unitful.NoUnits
    if i > 0
        println("   Rel. diff.: $diff")
    end
    println("Iteration $(i+1)...")
    diff > ϵ && i < maxiter
end
