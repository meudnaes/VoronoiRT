include("functions.jl")
include("characteristics.jl")

function J_λ(S_λ, α_tot, atmos; quadrature)

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

function Λ_iteration(ϵ, maxiter, atmos; quadrature)
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

    println(sum(isnan.(ε_λ[:, 2:end-1, 2:end-1])))

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
        J_new = J_λ(S_old, α_tot, atmos; quadrature=quadrature)
        J_test = J_new
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

function criterion(S_new, S_old, ϵ, i, maxiter)
    maximum(abs.((S_new[:,2:end-1,2:end-1] .- S_old[:,2:end-1,2:end-1])
                ./S_new[:,2:end-1,2:end-1])) > ϵ && i < maxiter
end
