# Push system closer to LTE, makes convergence easier.
# This is instead of developing operator splitting...
const BOOST = 2.0e9

"""
    calculate_C(atmos::Atmosphere,
                LTE_pops::Array{<:NumberDensity, 4})

Calculate collisional rates on the regular grid
"""
function calculate_C(atmos::Atmosphere,
                     LTE_pops::Array{<:NumberDensity, 4})

    nz, nx, ny, nl = size(LTE_pops)
    n_levels = nl - 1

    C = Array{Float64, 5}(undef, (n_levels+1, n_levels+1, nz, nx, ny))u"s^-1"

    # ionization
    for level = 1:n_levels
        C[level,n_levels+1,:,:,:] = Cij(level, n_levels+1, atmos.electron_density*1.0, atmos.temperature*1.0, LTE_pops)
        C[n_levels+1,level,:,:,:] = Cij(n_levels+1, level, atmos.electron_density*1.0, atmos.temperature*1.0, LTE_pops)
    end

    # bb transition
    # for l=1:n_levels-1
        # for u=(l+1):n_levels

    l = 1
    u = 2

    C[l,u,:,:,:] = Cij(l, u, atmos.electron_density*1.0, atmos.temperature*1.0, LTE_pops)
    C[u,l,:,:,:] = Cij(u, l, atmos.electron_density*1.0, atmos.temperature*1.0, LTE_pops)

        # end
    # end

    # Fill diagonal with zeros, because #undef does not like arithmetics
    for l=1:n_levels+1
        C[l,l,:,:,:] .= 0u"s^-1"
    end

    return C # ε ~ 0.1 at the upper parts, or on average, minimum should be 1e-2
end

"""
    calculate_C(sites::VoronoiSites,
                LTE_pops::Matrix{<:NumberDensity})

Calculate collisional rates on the regular grid
"""
function calculate_C(sites::VoronoiSites,
                     LTE_pops::Matrix{<:NumberDensity})

    n, nl = size(LTE_pops)
    n_levels = nl - 1

    C = Array{Float64, 3}(undef, (n_levels+1, n_levels+1, n))u"s^-1"

    # ionization
    for level = 1:n_levels
        C[level,n_levels+1,:] = Cij(level, n_levels+1, sites.electron_density*1.0, sites.temperature*1.0, LTE_pops)
        C[n_levels+1,level,:] = Cij(n_levels+1, level, sites.electron_density*1.0, sites.temperature*1.0, LTE_pops)
    end

    # bb transition
    # for l=1:n_levels-1
        # for u=(l+1):n_levels

    l = 1
    u = 2

    C[l,u,:] = Cij(l, u, sites.electron_density*1.0, sites.temperature*1.0, LTE_pops)
    C[u,l,:] = Cij(u, l, sites.electron_density*1.0, sites.temperature*1.0, LTE_pops)

        # end
    # end

    # Fill diagonal with zeros, because #undef does not like arithmetics
    for l=1:n_levels+1
        C[l,l,:] .= 0u"s^-1"
    end

    return C
end

"""
    calculate_R(atmos::Atmosphere,
                     line::HydrogenicLine,
                     J_λ::Array{<:UnitsIntensity_λ,4},
                     damping_λ::Array{<:Float64, 4},
                     LTE_pops::Array{<:NumberDensity, 4})

Calculcates the radiative rates on the regular grid
"""
function calculate_R(atmos::Atmosphere,
                     line::HydrogenicLine,
                     J_λ::Array{<:UnitsIntensity_λ,4},
                     damping_λ::Array{<:Float64, 4},
                     LTE_pops::Array{<:NumberDensity, 4})

    nz, nx, ny, nl = size(LTE_pops)
    n_levels = nl - 1

    R = Array{Float64, 5}(undef, (n_levels+1, n_levels+1, nz, nx, ny))u"s^-1"

    # ionization
    for level = 1:n_levels
        start = line.λidx[level+1]+1
        stop = line.λidx[level+2]
        σ = σic(level, line, line.λ[start:stop])
        G = Gij(level, n_levels+1, line.λ[start:stop], atmos.temperature*1.0, LTE_pops)

        R[level,n_levels+1,:,:,:] = Rij(J_λ[start:stop,:,:,:], σ, line.λ[start:stop])
        R[n_levels+1,level,:,:,:] = Rji(J_λ[start:stop,:,:,:], σ, G, line.λ[start:stop])
    end

    # bb transition
    # for l=1:n_levels-1
        # for u=(l+1):n_levels

    start = line.λidx[1]+1
    stop = line.λidx[2]

    l = 1
    u = 2

    σ = σij(l, u, line, line.λ[start:stop], damping_λ[start:stop,:,:,:])
    G = Gij(l, u, line.λ[start:stop], atmos.temperature*1.0, LTE_pops)

    R[l,u,:,:,:] = Rij(J_λ[start:stop,:,:,:], σ, line.λ[start:stop])
    R[u,l,:,:,:] = Rji(J_λ[start:stop,:,:,:], σ, G, line.λ[start:stop])

        # end
    # end

    # Fill diagonal with zeros, because #undef does not like arithmetics
    for l=1:n_levels+1
        R[l,l,:,:,:] .= 0u"s^-1"
    end

    return R
end

"""
    calculate_R(sites::VoronoiSites,
                line::HydrogenicLine,
                J_λ::Matrix{<:UnitsIntensity_λ},
                damping_λ::Matrix{<:Float64},
                LTE_pops::Matrix{<:NumberDensity})

Calculcates the radiative rates on the irregular grid
"""
function calculate_R(sites::VoronoiSites,
                     line::HydrogenicLine,
                     J_λ::Matrix{<:UnitsIntensity_λ},
                     damping_λ::Matrix{<:Float64},
                     LTE_pops::Matrix{<:NumberDensity})

    n, nl = size(LTE_pops)
    n_levels = nl - 1

    R = Array{Float64, 3}(undef, (n_levels+1, n_levels+1, n))u"s^-1"

    # ionization
    for level = 1:n_levels
        start = line.λidx[level+1]+1
        stop = line.λidx[level+2]
        σ = σic(level, line, line.λ[start:stop])
        G = Gij(level, n_levels+1, line.λ[start:stop], sites.temperature*1.0, LTE_pops)

        R[level,n_levels+1,:] = Rij(J_λ[start:stop,:], σ, line.λ[start:stop])
        R[n_levels+1,level,:] = Rji(J_λ[start:stop,:], σ, G, line.λ[start:stop])
    end

    # bb transition
    # for l=1:n_levels-1
        # for u=(l+1):n_levels

    start = line.λidx[1]+1
    stop = line.λidx[2]

    l = 1
    u = 2

    σ = σij(l, u, line, line.λ[start:stop], damping_λ[start:stop, :])
    G = Gij(l, u, line.λ[start:stop], sites.temperature*1.0, LTE_pops)

    R[l,u,:] = Rij(J_λ[start:stop,:], σ, line.λ[start:stop])
    R[u,l,:] = Rji(J_λ[start:stop,:], σ, G, line.λ[start:stop])

        # end
    # end

    # Fill diagonal with zeros, because #undef does not like arithmetics
    for l=1:n_levels+1
        R[l,l,:] .= 0u"s^-1"
    end

    return R
end


"""
    Rij(J::Array{<:UnitsIntensity_λ, 4},
        σij::Array{<:Unitful.Area, 4},
        λ::Vector{<:Unitful.Length})

Raditive rate for excitation transitions.
"""
function Rij(J::Array{<:UnitsIntensity_λ, 4},
             σij::Array{<:Unitful.Area, 4},
             λ::Vector{<:Unitful.Length})

    nλ, nz, nx, ny = size(J)
    R = Array{Float64, 3}(undef, (nz, nx, ny))u"s^-1"
    fill!(R,0.0u"s^-1")

    for l=1:(nλ-1)
        R += @. 2π/hc*((λ[l]*σij[l,:,:,:]*J[l,:,:,:] +
                     λ[l+1]*σij[l+1,:,:,:]*J[l+1,:,:,:])*(λ[l+1] - λ[l]))/1000
    end

    return R
end
function Rij(J::Matrix{<:UnitsIntensity_λ},
             σij::Matrix{<:Unitful.Area},
             λ::Vector{<:Unitful.Length})

    nλ, n = size(J)
    R = Vector{Float64}(undef, n)u"s^-1"
    fill!(R,0.0u"s^-1")

    for l=1:(nλ-1)
        R += @. 2π/hc*((λ[l]*σij[l,:]*J[l,:] +
                     λ[l+1]*σij[l+1,:]*J[l+1,:])*(λ[l+1] - λ[l]))/1000
    end

    return R
end

"""
    Rij(J::Array{<:UnitsIntensity_λ, 4},
        σij::Array{<:Unitful.Area, 1},
        λ::Vector{<:Unitful.Length})

Radiative rate for ionisation transitions.
"""
function Rij(J::Array{<:UnitsIntensity_λ, 4},
             σij::Vector{<:Unitful.Area},
             λ::Vector{<:Unitful.Length})

    nλ, nz, nx, ny = size(J)
    R = Array{Float64, 3}(undef, (nz, nx, ny))u"s^-1"
    fill!(R,0.0u"s^-1")

    for l=1:(nλ-1)
        R += @. 2π/hc*(λ[l]*σij[l]*J[l,:,:,:] +
                     λ[l+1]*σij[l+1]*J[l+1,:,:,:])*(λ[l+1] - λ[l])/1000
    end

    return R
end
function Rij(J::Matrix{<:UnitsIntensity_λ},
             σij::Vector{<:Unitful.Area},
             λ::Vector{<:Unitful.Length})

    nλ, n = size(J)
    R = Vector{Float64}(undef, n)u"s^-1"
    fill!(R,0.0u"s^-1")

    for l=1:(nλ-1)
        R += @. 2π/hc*(λ[l]*σij[l]*J[l,:] +
                     λ[l+1]*σij[l+1]*J[l+1,:])*(λ[l+1] - λ[l])/1000
    end

    return R
end

"""
    Rji(J::Array{<:UnitsIntensity_λ, 4},
        σij::Array{<:Unitful.Area, 4},
        Gij::Array{Float64, 4},
        λ::Vector{<:Unitful.Length})

Radiative rate for de-excitation transitions.
"""
function Rji(J::Array{<:UnitsIntensity_λ, 4},
             σij::Array{<:Unitful.Area, 4},
             Gij::Array{Float64, 4},
             λ::Vector{<:Unitful.Length})

    nλ, nz, nx, ny = size(J)
    R = Array{Float64, 3}(undef, (nz, nx, ny))u"s^-1"
    fill!(R,0.0u"s^-1")

    # Trapezoid rule
    for l=1:(nλ-1)
        R += @. 2π/hc*(σij[l,:,:,:]*Gij[l,:,:,:]*λ[l]*(2*h*c_0^2/λ[l]^5 + J[l,:,:,:] ) +
                    σij[l+1,:,:,:]*Gij[l+1,:,:,:]*λ[l+1]*(2*h*c_0^2 / λ[l+1]^5 + J[l+1,:,:,:] ))*(λ[l+1] - λ[l])
    end

    return R
end
function Rji(J::Matrix{<:UnitsIntensity_λ},
             σij::Matrix{<:Unitful.Area},
             Gij::Matrix{Float64},
             λ::Vector{<:Unitful.Length})

    nλ, n = size(J)
    R = Vector{Float64}(undef, n)u"s^-1"
    fill!(R,0.0u"s^-1")

    # Trapezoid rule
    for l=1:(nλ-1)
        R += @. 2π/hc*(σij[l,:]*Gij[l,:]*λ[l]*(2*h*c_0^2/λ[l]^5 + J[l,:] ) +
                    σij[l+1,:]*Gij[l+1,:]*λ[l+1]*(2*h*c_0^2 / λ[l+1]^5 + J[l+1,:] ))*(λ[l+1] - λ[l])
    end

    return R
end

"""
    Rji(J::Array{<:UnitsIntensity_λ, 4},
        σij::Array{<:Unitful.Area, 1},
        Gij::Array{Float64, 4},
        λ::Vector{<:Unitful.Length})

Radiative rate for recombination transitions.
"""
function Rji(J::Array{<:UnitsIntensity_λ, 4},
             σij::Vector{<:Unitful.Area},
             Gij::Array{Float64, 4},
             λ::Vector{<:Unitful.Length})

    nλ, nz, nx, ny = size(J)
    R = Array{Float64, 3}(undef, (nz, nx, ny))u"s^-1"
    fill!(R,0.0u"s^-1")

    # Trapezoid rule
    for l=1:(nλ-1)
        R += @. 2π/hc * (σij[l]*Gij[l,:,:,:]*λ[l]*(2*hc*c_0/λ[l]^5 + J[l,:,:,:]) +
                      σij[l+1]*Gij[l+1,:,:,:]*λ[l+1]*(2*hc*c_0/λ[l+1]^5 + J[l+1,:,:,:]))*(λ[l+1] - λ[l])
    end

    return R
end
function Rji(J::Matrix{<:UnitsIntensity_λ},
             σij::Vector{<:Unitful.Area},
             Gij::Matrix{Float64},
             λ::Vector{<:Unitful.Length})

    nλ, n = size(J)
    R = Vector{Float64}(undef, n)u"s^-1"
    fill!(R,0.0u"s^-1")

    # Trapezoid rule
    for l=1:(nλ-1)
        R += @. 2π/hc * (σij[l]*Gij[l,:]*λ[l]*(2*hc*c_0/λ[l]^5 + J[l,:]) +
                      σij[l+1]*Gij[l+1,:]*λ[l+1]*(2*hc*c_0/λ[l+1]^5 + J[l+1,:]))*(λ[l+1] - λ[l])
    end

    return R
end

"""
    σij(i::Integer,
        j::Integer,
        atom::Atom,
        λ::Vector{<:Unitful.Length})

Calculates the bound-bound cross section.
"""
function σij(i::Integer,
             j::Integer,
             line::HydrogenicLine,
             λ::Vector{<:Unitful.Length},
             damping_λ::Array{<:Float64, 4})

    λ0 = line.λ0
    σ_constant = h*c_0/(4π*λ0) * line.Bij
    nλ = length(λ)
    nz, nx, ny = size(line.ΔD)
    σ = Array{Unitful.Area, 4}(undef, nλ, nz, nx, ny)

    for l=1:nλ
        v = (λ[l] - λ0) ./ line.ΔD
        profile = voigt_profile.(damping_λ[l, :, :, :], ustrip(v), line.ΔD)
        σ[l,:,:,:] = σ_constant .* profile
    end

    return σ
end
function σij(i::Integer,
             j::Integer,
             line::HydrogenicLine,
             λ::Vector{<:Unitful.Length},
             damping_λ::Matrix{<:Float64})

    λ0 = line.λ0
    σ_constant = h*c_0/(4π*λ0) * line.Bij
    nλ = length(λ)
    n = length(line.ΔD)
    σ = Matrix{Float64}(undef, (nλ, n))u"m^2"

    for l=1:nλ
        v = (λ[l] - λ0) ./ line.ΔD
        profile = voigt_profile.(damping_λ[l, :], ustrip(v), line.ΔD)
        σ[l, :] = σ_constant .* profile
    end

    return σ
end

"""
    σic(i::Integer,
        atom::Atom,
        λ::Vector{<:Unitful.Length})

Calculates the bound-free cross-section.
"""
function σic(i::Integer,
             line::HydrogenicLine,
             λ::Vector{<:Unitful.Length},
             λ_edge=nothing)

    if λ_edge == nothing
        λ_edge = λ[end]
    end
    λ3_ratio = (λ ./ λ_edge).^3
    n_eff = sqrt(E_∞ / (line.χj - line.χi)) |>u"J/J" # should be χu - χl
    charge = line.Z
    σ_constant = (4 * e^2 / (3 * π * sqrt(3) * ε_0 * m_e * c_0^2 * R_∞)) |> u"m^2"

    σ = (σ_constant * charge^4 * n_eff * λ3_ratio .* gaunt_bf.(λ, charge, n_eff))

    return σ
end

"""
    Gij(i::Integer,
        j::Integer,
        λ::Vector{<:Unitful.Length},
        temperature::Array{<:Unitful.Temperature, 3},
        LTE_populations::Array{<:NumberDensity, 4})

Factor to go into the de-excitation and recombination expressions.
"""
function Gij(i::Integer,
             j::Integer,
             λ::Vector{<:Unitful.Length},
             temperature::Array{<:Unitful.Temperature, 3},
             LTE_populations::Array{<:NumberDensity, 4})

    nλ = length(λ)
    nz, nx, ny = size(temperature)
    G = Array{Float64, 4}(undef, nλ, nz, nx, ny)

    n_ratio = LTE_populations[:,:,:,i] ./LTE_populations[:,:,:,j]

    for l=1:nλ
        G[l,:,:,:] = @. n_ratio * exp(- hc / (k_B*λ[l]*temperature))
    end

    return G
end
function Gij(i::Integer,
             j::Integer,
             λ::Vector{<:Unitful.Length},
             temperature::Vector{<:Unitful.Temperature},
             LTE_populations::Matrix{<:NumberDensity})

    nλ = length(λ)
    n = length(temperature)
    G = Matrix{Float64}(undef, (nλ, n))

    n_ratio = LTE_populations[:,i] ./LTE_populations[:,j]

    for l=1:nλ
        G[l,:] = @. n_ratio * exp(- hc / (k_B*λ[l]*temperature))
    end

    return G
end

"""
    Cij(i::Integer,
        j::Integer,
        electron_density::Array{<:NumberDensity,3},
        temperature::Array{<:Unitful.Temperature,3},
        LTE_populations::Array{<:NumberDensity,4})

Calculates the collisional rates for all possible
two-level atom transitions.
"""
function Cij(i::Integer,
             j::Integer,
             electron_density::Array{<:NumberDensity,3},
             temperature::Array{<:Unitful.Temperature,3},
             LTE_populations::Array{<:NumberDensity,4})

    ionisation_level = size(LTE_populations)[end]

    # If UP
    if i < j
        if j < ionisation_level
            C = coll_exc_hydrogen_johnson.(i, j, electron_density, temperature)
        elseif j == ionisation_level
            C = coll_ion_hydrogen_johnson.(i, electron_density, temperature)
        end

    # If DOWN
    elseif i > j
        if i < ionisation_level
            C = coll_exc_hydrogen_johnson.(j, i, electron_density, temperature)
        elseif i == ionisation_level
            C = coll_ion_hydrogen_johnson.(j, electron_density, temperature)
        end
        C = @. C * ( LTE_populations[:,:,:,j] / LTE_populations[:,:,:,i] )
    end

    return C.*BOOST
end
function Cij(i::Integer,
             j::Integer,
             electron_density::Array{<:NumberDensity,1},
             temperature::Array{<:Unitful.Temperature,1},
             LTE_populations::Array{<:NumberDensity,2})

    ionisation_level = size(LTE_populations)[end]

    # If UP
    if i < j
        if j < ionisation_level
            C = coll_exc_hydrogen_johnson.(i, j, electron_density, temperature)
        elseif j == ionisation_level
            C = coll_ion_hydrogen_johnson.(i, electron_density, temperature)
        end

    # If DOWN
    elseif i > j
        if i < ionisation_level
            C = coll_exc_hydrogen_johnson.(j, i, electron_density, temperature)
        elseif i == ionisation_level
            C = coll_ion_hydrogen_johnson.(j, electron_density, temperature)
        end
        C = @. C * ( LTE_populations[:,j] / LTE_populations[:,i] )
    end

    return C.*BOOST
end


"""
    gaunt_bf(charge::Int, n_eff::Number, λ::Unitful.Length)::Float64

Compute bound-free Gaunt factor for a given charge, effective principal
quantum number and wavelength λ. Taken from RH. Formula from
[Seaton (1960), Rep. Prog. Phys. 23, 313](https://ui.adsabs.harvard.edu/abs/1960RPPh...23..313S/abstract),
page 316. Recipes from Seaton. (This is copied from github.com/tiagopereira/Transparency.jl)
"""
function gaunt_bf(λ::Unitful.Length,
                  charge::Real,
                  n_eff::Real)::Float64

    x = ustrip(1 / (λ * R_∞ * charge^2) |> u"m/m")
    x3 = x^(1/3)
    nsqx = 1 / (n_eff^2 * x)
    g_bf = 1 + 0.1728 * x3 * (1 - 2 * nsqx) - 0.0496 * x3^2 * (1 - (1 - nsqx) * 0.66666667 * nsqx)
    @assert g_bf >= 0 "gaunt_bf negative, calculation will not be reliable"
    return g_bf
end
