"""
    zero_radiation_populations(line::HydrogenicLine,
                               atmos::Atmosphere)

For a given atom density, calculate the populations according to zero-radiation,
on a regular grid.
"""
function zero_radiation_populations(line::HydrogenicLine,
                                    atmos::Atmosphere)

    nz, nx, ny = size(atmos.temperature)
    nλ = length(line.λ)

    J_zero = zeros(Float64, (nλ, nz, nx, ny))u"kW*m^-2*nm^-1"

    LTE_pops = LTE_populations(line, atmos)

    γ = γ_constant(line,
                   atmos.temperature,
                   (LTE_pops[:, :, :, 1] .+ LTE_pops[:, :, :, 2]),
                   atmos.electron_density)

    damping_λ = Array{Float64, 4}(undef, size(J_zero))
    for l in eachindex(line.λ)
       damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    R_zero = calculate_R(atmos, line, J_zero, damping_λ, LTE_pops)
    C_zero = calculate_C(atmos, LTE_pops)

    populations = get_revised_populations(R_zero, C_zero, atmos.hydrogen_populations)

    return populations
end

"""
    zero_radiation_populations(line::HydrogenicLine,
                               sites::VoronoiSites)

For a given atom density, calculate the populations according to zero-radiation,
on an irregular grid.
"""
function zero_radiation_populations(line::HydrogenicLine,
                                    sites::VoronoiSites)

    nλ = length(line.λ)

    J_zero = zeros(Float64, (nλ, sites.n))u"kW*m^-2*nm^-1"

    LTE_pops = LTE_populations(line, sites)

    γ = γ_constant(line,
                   sites.temperature,
                   (LTE_pops[:, 1] .+ LTE_pops[:, 2]),
                   sites.electron_density)

    damping_λ = Matrix{Float64}(undef, size(J_zero))
    for l in eachindex(line.λ)
       damping_λ[l, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    R_zero = calculate_R(sites, line, J_zero, damping_λ, LTE_pops)
    C_zero = calculate_C(sites, LTE_pops)

    populations = get_revised_populations(R_zero, C_zero, sites.hydrogen_populations)

    return populations
end

"""
    LTE_populations(line::HydrogenicLine,
                    atmos::Atmosphere)

Given the atom density, calculate the atom populations according to LTE, on a
regular grid.
"""
function LTE_populations(line::HydrogenicLine,
                         atmos::Atmosphere)

    χ = [line.χi, line.χj, line.χ∞]
    # Ionised hydrogen -> g = 1
    g = [line.gi, line.gj, 1]
    atom_density = atmos.hydrogen_populations*1.0
    nz, nx, ny = size(atom_density)

    n_levels = 3
    n_relative = ones(Float64, nz, nx, ny, n_levels)

    saha_const = (k_B / h) * (2π * m_e) / h
    saha_factor = 2 * ((saha_const * atmos.temperature*1.0).^(3/2) ./ atmos.electron_density*1.0) .|> u"m/m"

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

"""
    LTE_populations(line::HydrogenicLine,
                    sites::VoronoiSites)

Given the atom density, calculate the atom populations according to LTE, on an
irregular grid.
"""
function LTE_populations(line::HydrogenicLine,
                         sites::VoronoiSites)

    χ = [line.χi, line.χj, line.χ∞]
    # Ionised hydrogen -> g = 1
    g = [line.gi, line.gj, 1]
    atom_density = sites.hydrogen_populations
    n = length(atom_density)

    n_levels = 3
    n_relative = ones(Float64, n, n_levels)

    saha_const = (k_B / h) * (2π * m_e) / h
    saha_factor = 2 * ((saha_const * sites.temperature).^(3/2) ./ sites.electron_density) .|> u"m/m"

    for i=2:n_levels
        ΔE = χ[i] - χ[1]
        n_relative[:, i] = g[i] / g[1] * exp.(-ΔE ./ (k_B * sites.temperature))
    end

    # Last level is ionised stage (H II)
    n_relative[:,n_levels] .*= saha_factor
    n_relative[:,1] = 1 ./ sum(n_relative, dims=2)[:,1]
    n_relative[:,2:end] .*= n_relative[:,1]

    return n_relative .* atom_density
end

"""
    get_revised_populations(R::Array{<:Unitful.Frequency, 5},
                            C::Array{<:Unitful.Frequency, 5},
                            atom_density::Array{<:NumberDensity,3})

Update populations on a regular grid with Statistical Equilibrium
"""
function get_revised_populations(R::Array{<:Unitful.Frequency, 5},
                                 C::Array{<:Unitful.Frequency, 5},
                                 atom_density::Array{<:NumberDensity,3})

    # Total rates
    P = R .+ C

    n_levels = size(P)[1] - 1
    nz, nx, ny = size(atom_density)

    A = Array{Float64, 5}(undef, (n_levels, n_levels, nz, nx, ny))u"s^-1"
    b = Array{Float64, 4}(undef, (n_levels, nz, nx, ny))u"s^-1*m^-3"
    populations = Array{Float64, 4}(undef, (nz, nx, ny, n_levels + 1))u"m^-3"

    for r = 1:n_levels
        A[r, r, :, :, :] = P[1, r+1, :, :, :] .+ P[r+1, 1, :, :, :]
        for c = setdiff(1:n_levels, r)
            A[r, c, :, :, :] = P[1, r+1, :, :, :] .- P[c+1, r+1, :, :, :]
            A[r, r, :, :, :] .+= P[r+1, c+1, :, :, :]
        end

        b[r, :, :, :] = atom_density .* P[1, r+1, :, :, :]
    end

    for k = 1:nz
        for i = 1:nx
            for j = 1:ny
                populations[k, i, j, 2:end] .= inv(A[:, :, k, i, j])*b[:, k, i, j]
            end
        end
    end

    populations[: ,: ,: , 1] = atom_density .- sum(populations[:, :, :, 2:end], dims=4)[:, :, :, 1]

    return populations
end

"""
    get_revised_populations(R::Array{<:Unitful.Frequency, 3},
                            C::Array{<:Unitful.Frequency, 3},
                            atom_density::Vector{<:NumberDensity})

Update populations on an irregular grid with Statistical Equilibrium
"""
function get_revised_populations(R::Array{<:Unitful.Frequency, 3},
                                 C::Array{<:Unitful.Frequency, 3},
                                 atom_density::Vector{<:NumberDensity})

    P = R .+ C

    n_levels = size(P)[1] - 1
    n = length(atom_density)

    A = Array{Float64, 3}(undef, (n_levels, n_levels, n))u"s^-1"
    b = Matrix{Float64}(undef, (n_levels, n))u"s^-1*m^-3"
    populations = Matrix{Float64}(undef, (n, n_levels + 1))u"m^-3"

    for r=1:n_levels
        A[r,r,:] = P[1,r+1,:] .+ P[r+1,1,:]
        for c=setdiff(1:n_levels, r)
            A[r,c,:] = P[1,r+1,:] .- P[c+1,r+1,:]
            A[r,r,:] .+= P[r+1,c+1,:]
        end

        b[r,:,:,:] = atom_density .* P[1,r+1,:,:,:]
    end

    for ii=1:n
        populations[ii, 2:end] .= inv(A[:, :, ii])*b[:, ii]
    end

    populations[:, 1] = atom_density .- sum(populations[:, 2:end], dims=2)[:,1]

    return populations
end
