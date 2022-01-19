include("line.jl")
include("functions.jl")
include("voronoi_utils.jl")
include("characteristics.jl")
include("irregular_ray_tracing.jl")

function J_λ_regular(S_λ::AbstractArray, α_tot::AbstractArray, atmos::Atmosphere, quadrature::String)

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

function J_λ_voronoi(S_λ::AbstractArray, α_tot::AbstractArray, sites::VoronoiSites, quadrature::String)

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

function Λ_regular(ϵ::AbstractFloat, maxiter::Integer, atmos::Atmosphere, quadrature::String)
    # choose a wavelength
    λ = 500u"nm"  # nm

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_tot = α_cont.(λ, atmos.temperature*1.0, atmos.electron_density*1.0,
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

function Λ_regular(ϵ::AbstractFloat, maxiter::Integer, atmos::Atmosphere, line::HydrogenicLine, quadrature::String)
    # choose a wavelength
    λ = 500u"nm"  # nm

    LTE_pops = LTE_populations(line, atmos)

    h_ground_density = LTE_pops[:, :, :, 1] .+ LTE_pops[:, :, :, 2]
    proton_density = LTE_pops[:, :, :, 3]

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_tot = α_cont.(λ, atmos.temperature*1.0, atmos.electron_density*1.0,
                    h_ground_density*1.0, proton_density*1.0)

    α_a = α_absorption.(λ, atmos.temperature*1.0, atmos.electron_density*1.0,
                        h_ground_density*1.0, proton_density*1.0)


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

function Λ_voronoi(ϵ::AbstractFloat, maxiter::Integer, sites::VoronoiSites, quadrature::String)
    println("---Iterating---")
    # choose a wavelength
    λ = 500u"nm"  # nm

    # Only continuum
    η_ν = 0

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_tot = α_cont.(λ, sites.temperature*1.0, sites.electron_density*1.0,
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

function criterion(S_new, S_old, ϵ, i, maxiter, indcs)
    diff = maximum(abs.((S_new[indcs] .- S_old[indcs])./S_new[indcs])) |> Unitful.NoUnits
    if i > 0
        println(diff)
    end
    diff > ϵ && i < maxiter
end
