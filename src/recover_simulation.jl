include("io.jl")
include("line.jl")
include("plot_utils.jl")
include("functions.jl")
include("atmosphere.jl")
include("lambda_iteration.jl")

global nλ_bb = 50
global nλ_bf = 20

function recover_regular(ϵ::AbstractFloat,
                         maxiter::Integer,
                         quadrature::String,
                         DATA::String)

    println("---Recovering simulation---")
    atmos, S_new, populations = read_quantities(DATA; periodic=true)
    line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., atmos.temperature)

    # LTE populations
    LTE_pops = LTE_populations(line, atmos)

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

    S_old = zero(B_0)

    C = calculate_C(atmos, LTE_pops)

    local convergence
    h5open(DATA, "r") do file
        convergence = read(file, "convergence")[:]
    end

    local i
    for j in 1:length(convergence)
        if convergence[j] == 0
            i = j
            break
        end
    end

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

    h5open(DATA, "r+") do file
       file["time"][:] = 0.8888888888
    end

    return J_new, S_new, α_cont, populations
end

function recover_voronoi(ϵ::AbstractFloat,
                         maxiter::Integer,
                         quadrature::String,
                         DATA::String)

    println("---Recovering simulation---")
    sites, S_new, populations = read_sites(DATA)
    line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., sites.temperature)

    println(size(S_new))
    println(length(line.λ))

    return

    # Start in LTE
    LTE_pops = LTE_populations(line, sites)

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

    S_old = zero(B_0)

    # destruction probability (Should I include line???)
    ελ = destruction(LTE_pops, sites.electron_density, sites.temperature, line)
    println("Minimum $(minimum(ελ)) destruction probability")

    C = calculate_C(sites, LTE_pops)

    local convergence
    h5open(DATA, "r") do file
        convergence = read(file, "convergence")[:]
    end

    local i
    for j in 1:length(convergence)
        if convergence[j] == 0
            i = j
            break
        end
    end

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

    h5open(DATA, "r+") do file
       file["time"][:] = 0.8888888888
    end

    return J_new, S_new, α_cont, populations
end

"""
    read_sites(DATA::String)

read quantities from irregular grid simulation from a hdf5 file
"""
function read_sites(DATA::String)
    local positions, temperature, electron_density, hydrogen_populations, S_λ, populations
    local velocity_z, velocity_x, velocity_y, boundaries
    h5open(DATA, "r") do file
        positions = read(file, "positions")[:, :]*u"m"
        boundaries = read(file, "boundaries")[:]u"m"

        temperature = read(file, "temperature")[:]*u"K"

        electron_density = read(file, "electron_density")[:]*u"m^-3"
        hydrogen_populations = read(file, "hydrogen_populations")[:]*u"m^-3"

        velocity_z = read(file, "velocity_z")[:]u"m*s^-1"
        velocity_x = read(file, "velocity_x")[:]u"m*s^-1"
        velocity_y = read(file, "velocity_y")[:]u"m*s^-1"

        S_λ = read(file, "source_function")[:, :]*u"kW*m^-2*nm^-1"
        populations = read(file, "populations")[:, :]*u"m^-3"
    end

    sites_file = "../data/sites_recover.txt"
    neighbours_file = "../data/neighbours_recover.txt"
    # write sites to file
    write_arrays(positions[2, :],
                 positions[3, :],
                 positions[1, :],
                 sites_file)

    x_min = ustrip(boundaries[3])
    x_max = ustrip(boundaries[4])
    y_min = ustrip(boundaries[5])
    y_max = ustrip(boundaries[6])
    z_min = ustrip(boundaries[1])
    z_max = ustrip(boundaries[2])

    # export sites to voro++, and compute grid information
    println("---Preprocessing grid---")


    # compute neigbours
    run(`./voro.sh $sites_file $neighbours_file
                   $(x_min) $(x_max)
                   $(y_min) $(y_max)
                   $(z_min) $(z_max)`)

    n_sites = length(positions[1,:])

    # Voronoi grid
    sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions,
                                   x_min*1u"m", x_max*1u"m",
                                   y_min*1u"m", y_max*1u"m")...,
                         temperature,
                         electron_density,
                         hydrogen_populations,
                         velocity_z,
                         velocity_x,
                         velocity_y,
                         z_min*1u"m", z_max*1u"m",
                         x_min*1u"m", x_max*1u"m",
                         y_min*1u"m", y_max*1u"m",
                         n_sites)


    return sites, S_λ, populations
end

recover_voronoi(1e-3, 100, "../quadratures/ul7n12.dat",
                "../data/voronoi_ul7n12_extinction_3e5_recover.h5")
