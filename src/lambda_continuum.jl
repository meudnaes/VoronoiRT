include("io.jl")
include("functions.jl")
include("radiation.jl")
include("voronoi_utils.jl")
include("characteristics.jl")
include("irregular_ray_tracing.jl")


function J_λ_regular(S_λ::AbstractArray,
                     α_cont::AbstractArray,
                     atmos::Atmosphere,
                     quadrature::String)

    # Ω = (θ, φ), space angle
    weights, θ_array, ϕ_array, n_points = read_quadrature(quadrature)

    J = zero(S_λ)

    for i in 1:n_points
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        if θ > 90
            I_0 = blackbody_λ.(500u"nm", atmos.temperature[1,:,:])
            J += weights[i]*short_characteristics_up(k, S_λ, α_cont, atmos, I_0=I_0)
        elseif θ < 90
            I_0 = zero(S_λ[1, :, :])
            J += weights[i]*short_characteristics_down(k, S_λ, α_cont, atmos, I_0=I_0)
        end
    end
    return J
end


function J_λ_voronoi(S_λ::AbstractArray,
                     α_cont::AbstractArray,
                     sites::VoronoiSites,
                     quadrature::String)

    # Ω = (θ, φ), space angle
    weights, θ_array, ϕ_array, n_points = read_quadrature(quadrature)

    J = zero(S_λ)

    n_sweeps = 3
    Nran = 1

    for i in 1:n_points
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        if θ > 90
            bottom_layer = sites.layers_up[2] - 1
            bottom_layer_idx = sites.perm_up[1:bottom_layer]
            I_0 = blackbody_λ.(500u"nm", sites.temperature[bottom_layer_idx])
            J += weights[i]*Delaunay_upII(k, S_λ, I_0, α_cont, sites, n_sweeps)
        elseif θ < 90
            top_layer = sites.layers_down[2] - 1
            I_0 = zeros(top_layer)u"kW*nm^-1*m^-2"
            J += weights[i]*Delaunay_downII(k, S_λ, I_0, α_cont, sites, n_sweeps)
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

    # Lte populations
    LTE_pops = LTE_populations(atmos)

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_s = α_scattering(λ, atmos.temperature, atmos.electron_density*1.0,
                       LTE_pops[:,:,:,1])

    α_a = α_absorption.(λ, atmos.temperature, atmos.electron_density*1.0,
                        LTE_pops[:,:,:,1].+LTE_pops[:,:,:,2],
                        LTE_pops[:,:,:,3])

    α_cont = α_s + α_a
    # destruction
    ε_λ = α_a./α_cont

    thick = ε_λ .> 1e-4

    # Start with the source function as the Planck function
    B_0 = blackbody_λ.(λ, atmos.temperature)
    S_new = B_0

    S_old = zero(S_new)

    i=0

    local J_new
    # check where ε < 1e-2, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter, thick)
        S_old = copy(S_new)
        J_new = J_λ_regular(S_old, α_cont, atmos, quadrature)
        S_new = (1 .- ε_λ).*J_new .+ ε_λ.*B_0
        i+=1
    end

    if i == maxiter
        println("Did not converge inside scope")
        return J_new, S_new, α_cont
    end

    println("Converged in $i iterations")

    return J_new, S_new, α_cont
end

function Λ_voronoi(ϵ::AbstractFloat,
                   maxiter::Integer,
                   sites::VoronoiSites,
                   quadrature::String)
    println("---Iterating---")

    # choose a wavelength
    λ = 500u"nm"  # nm

    # Lte populations
    LTE_pops = LTE_populations(sites)

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_s = α_scattering(λ, sites.temperature, sites.electron_density*1.0,
                       LTE_pops[:,1])

    α_a = α_absorption.(λ, sites.temperature, sites.electron_density*1.0,
                        LTE_pops[:,1].+LTE_pops[:,2],
                        LTE_pops[:,3])

    α_cont = α_s + α_a

    # destruction
    ε_λ = α_a ./ α_cont

    thick = ε_λ .> 1e-4

    # Start with the source function as the Planck function
    B_0 = blackbody_λ.(λ, sites.temperature)
    S_new = B_0

    S_old = zero(S_new)

    i=0

    local J_new
    # check where ε < 1e-2, cut above heights
    while criterion(S_new, S_old, ϵ, i, maxiter, thick)
        S_old = copy(S_new)
        J_new = J_λ_voronoi(S_old, α_cont, sites, quadrature)
        S_new = (1 .- ε_λ).*J_new .+ ε_λ.*B_0
        i+=1
    end

    if i == maxiter
        println("Did not converge inside scope")
        return J_new, S_new, α_cont
    end

    println("Converged in $i iterations")

    return J_new, S_new, α_cont
end

function criterion(S_new::Array{<:UnitsIntensity_λ, 3},
                   S_old::Array{<:UnitsIntensity_λ, 3},
                   ϵ::Float64,
                   i::Int,
                   maxiter::Int,
                   indcs)

    diff = maximum(abs.(1 .- S_old[indcs]./S_new[indcs])) |> Unitful.NoUnits
    if isnan(diff)
        println("NaN DIFF!, index $l")
    end
    if i > 0
        println("   Rel. diff.: $diff")
    end
    println("Iteration $(i+1)...")

    diff > ϵ && i < maxiter
end

function criterion(S_new::Vector{<:UnitsIntensity_λ},
                   S_old::Vector{<:UnitsIntensity_λ},
                   ϵ::Float64,
                   i::Int,
                   maxiter::Int,
                   indcs)

    diff = maximum(abs.(1 .- S_old[indcs]./S_new[indcs])) |> Unitful.NoUnits
    if isnan(diff)
        println("NaN DIFF!, index $l")
    end
    if i > 0
        println("   Rel. diff.: $diff")
    end
    println("Iteration $(i+1)...")

    diff > ϵ && i < maxiter
end
