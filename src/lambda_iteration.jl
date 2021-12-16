include("functions.jl")
include("voronoi_utils.jl")
include("characteristics.jl")
include("irregular_ray_tracing.jl")

function J_λ_regular(S_λ::AbstractArray, α_tot::AbstractArray, atmos::Atmosphere; quadrature)

    # Ω = (θ, φ), space angle
    weights, θ_array, ϕ_array, n_points = read_quadrature(quadrature)

    J = zero(S_λ)

    for i in 1:n_points
        if θ_array[i] > 90
            J += weights[i]*short_characteristics_up(θ_array[i], ϕ_array[i], S_λ,
                                                     α_tot, atmos, degrees=true)
        elseif θ_array[i] < 90
            J += weights[i]*short_characteristics_down(θ_array[i], ϕ_array[i], S_λ,
                                                       α_tot, atmos, degrees=true)
        end
    end
    return J
end

function J_λ_voronoi(S_λ::AbstractArray, α_tot::AbstractArray, sites::VoronoiSites; quadrature)

    # Ω = (θ, φ), space angle
    weights, θ_array, ϕ_array, n_points = read_quadrature(quadrature)

    J = zero(S_λ)

    for i in 1:n_points
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        k = -[cos(θ), cos(ϕ)*sin(θ), sin(ϕ)*sin(θ)]
        if θ > 90
            J += weights[i]*SC_Delaunay_up(sites, I_0, S_0, α_0,
                                           S, α, k, n_sweeps, Nran)
        elseif θ < 90
            J += weights[i]*SC_Delaunay_down(sites, I_0, S_0, α_0,
                                             S, α, k, n_sweeps, Nran)
        end
    end
    return J
end

function Λ_regular(ϵ::AbstractFloat, maxiter::Integer, atmos::Atmosphere, quadrature)
    # choose a wavelength
    λ = 500u"nm"  # nm

    # Only continuum
    η_ν = 0

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_tot = α_cont.(λ, atmos.temperature, atmos.electron_density,
                atmos.hydrogen_populations, atmos.hydrogen_populations)

    α_a = α_absorption.(λ, atmos.temperature, atmos.electron_density,
                     atmos.hydrogen_populations, atmos.hydrogen_populations)

    # destruction
    ε_λ = α_a ./ α_tot

    # Start with the source function as the Planck function
    S_new = B_λ.(λ, atmos.temperature)
    B_0 = S_new

    S_old = zero(S_new)

    i=0

    local J_new
    # check where ε < 1e-2, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter)
        print("Iteration $(i+1)\r")
        S_old = S_new
        J_new = J_λ(S_old, α_tot, atmos, quadrature)
        J_test = J_new
        S_new = (1 .- ε_λ).*J_new .+ ε_λ.*B_0
        i+=1
    end

    if i == maxiter
        println("Did not converge inside scope")
        return J_new, S_new
    end

    println("Converged in $i iterations")

    return J_new, S_new, α_tot
end

function Λ_voronoi(ϵ::AbstractFloat, maxiter::Integer, sites::VoronoiSites, quadrature)
    # choose a wavelength
    λ = 500u"nm"  # nm

    # Only continuum
    η_ν = 0

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_tot = α_cont.(λ, atmos.temperature, atmos.electron_density,
                atmos.hydrogen_populations, atmos.hydrogen_populations)

    α_a = α_absorption.(λ, atmos.temperature, atmos.electron_density,
                     atmos.hydrogen_populations, atmos.hydrogen_populations)

    # destruction
    ε_λ = α_a ./ α_tot

    # Start with the source function as the Planck function
    S_new = B_λ.(λ, atmos.temperature)
    B_0 = S_new

    S_old = zero(S_new)

    i=0

    local J_new
    # check where ε < 1e-2, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter)
        print("Iteration $(i+1)\r")
        S_old = S_new
        J_new = J_λ(S_old, α_tot, sites, quadrature)
        J_test = J_new
        S_new = (1 .- ε_λ).*J_new .+ ε_λ.*B_0
        i+=1
    end

    if i == maxiter
        println("Did not converge inside scope")
        return J_new, S_new
    end

    println("Converged in $i iterations")

    return J_new, S_new, α_tot
end

function criterion(S_new, S_old, ϵ, i, maxiter)
    maximum(abs.((S_new[:,2:end-1,2:end-1] .- S_old[:,2:end-1,2:end-1])
                ./S_new[:,2:end-1,2:end-1])) > ϵ && i < maxiter
end
