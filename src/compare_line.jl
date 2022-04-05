using Plots

include("io.jl")
# include("line.jl")
include("plot_utils.jl")
include("functions.jl")
include("atmosphere.jl")
include("lambda_iteration.jl")

global my_seed = 2022
Random.seed!(my_seed)

function sample_from_destruction(atmos::Atmosphere)

    nλ_bb = 0
    nλ_bf = 0
    line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., atmos.temperature)
    LTE_pops = LTE_populations(line, atmos)
    ελ = destruction(LTE_pops, atmos.electron_density, atmos.temperature, line)

    return ustrip.(ελ)
end

function compare(DATA, quadrature)
    maxiter = 100
    println("---Iterating maximum $maxiter iterations---")
    println("--- ! Boosting collisional rates ! ---")
    ϵ = 1e-3

    n_skip = 1

    nλ_bb = 50
    nλ_bf = 20

    function regular()

        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=4)...)
        line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., atmos.temperature)

        nx = length(atmos.x[2:end-1])
        ny = length(atmos.y[2:end-1])
        nz = length(atmos.z[2:end-1])

        println("sites: $(nx*ny*nz)")

        REGULAR_DATA = "../data/test_LTE_pops.h5"
        # "regular_ul7n12_half_C_2e9_single.h5"

        create_output_file(REGULAR_DATA, length(line.λ), size(atmos.temperature[:, 2:end-1, 2:end-1]), maxiter)
        write_to_file(atmos, REGULAR_DATA, ghost_cells=true)
        write_to_file(nλ_bb, "n_bb", REGULAR_DATA)
        write_to_file(nλ_bf, "n_bf", REGULAR_DATA)
        write_to_file(line, RTEGULAR_DATA)

        (J_mean, S_λ, α_cont, populations), time = @timed Λ_regular(ϵ, maxiter, atmos, line, quadrature, REGULAR_DATA)

        h5open(REGULAR_DATA, "r+") do file
           file["time"][:] = time
        end

        return
    end


    function voronoi()

        atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=n_skip)...)

        nx = length(atmos.x)
        nz = length(atmos.z)
        ny = length(atmos.y)

        n_sites = 3_000_000
        # 286720# floor(Int, nz*nx*ny)

        positions = sample_from_total_extinction(atmos, n_sites)
        # rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))
        # sample_from_destruction(atmos, n_sites)
        # sample_from_extinction(atmos, 121.562u"nm", n_sites)

        sites_file = "../data/sites_compare.txt"
        neighbours_file = "../data/neighbours_compare.txt"
        # write sites to file
        write_arrays(positions[2, :],
                     positions[3, :],
                     positions[1, :],
                     sites_file)

        x_min = ustrip(atmos.x[1])
        x_max = ustrip(atmos.x[end])
        y_min = ustrip(atmos.y[1])
        y_max = ustrip(atmos.y[end])
        z_min = ustrip(atmos.z[1])
        z_max = ustrip(atmos.z[end])

        # export sites to voro++, and compute grid information
        println("---Preprocessing grid---")


        # compute neigbours
        run(`./voro.sh $sites_file $neighbours_file
                       $(x_min) $(x_max)
                       $(y_min) $(y_max)
                       $(z_min) $(z_max)`)

        # Voronoi grid
        sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions)...,
                             _initialise(positions, atmos)...,
                             z_min*1u"m", z_max*1u"m",
                             x_min*1u"m", x_max*1u"m",
                             y_min*1u"m", y_max*1u"m",
                             n_sites)

        println("--grid calculations done---")
        run(`rm ../data/$sites_file ../data/$neighbours_file`)

        line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., sites.temperature)

        VORONOI_DATA = "../data/total_ext_3e6_copy.h5"

        create_output_file(VORONOI_DATA, length(line.λ), n_sites, maxiter)
        write_to_file(nλ_bb, "n_bb", VORONOI_DATA)
        write_to_file(nλ_bf, "n_bf", VORONOI_DATA)
        write_to_file(sites, VORONOI_DATA)
        write_to_file(line, VORONOI_DATA)

        (J_mean, S_λ, α_cont, populations), time = @timed Λ_voronoi(ϵ, maxiter, sites, line, quadrature, VORONOI_DATA)

        h5open(VORONOI_DATA, "r+") do file
           file["time"][:] = time
        end

        return
    end

    # regular();
    # voronoi();
end

function LTE_line(DATA)
    θ = 180.0
    ϕ = 0.0

    nλ_bb = 50
    nλ_bf = 20

    atmos = Atmosphere(get_atmos(DATA; periodic=true)...)
    line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., atmos.temperature)

    # LTE populations
    LTE_pops = LTE_populations(line, atmos)
    populations = copy(LTE_pops)

    # LTE_data = "../data/LTE_data.h5"
    # create_output_file(LTE_data, length(line.λ),  size(atmos.temperature[:, 2:end-1, 2:end-1]), 1)
    # write_to_file(nλ_bb, "n_bb", LTE_data)
    # write_to_file(nλ_bf, "n_bf", LTE_data)
    # write_to_file(atmos, LTE_data, ghost_cells=true)
    # write_to_file(LTE_pops[:, 2:end-1, 2:end-1, :], LTE_data)

    # The source function is the Planck function
    B_0 = Array{Float64, 4}(undef, (length(line.λ), size(atmos.temperature)...))u"kW*m^-2*nm^-1"
    for l in eachindex(line.λ)
        B_0[l, :, :, :] = B_λ.(line.λ[l], atmos.temperature*1.0)
    end
    S_λ = copy(B_0)

    # Find continuum extinction and absorption extinction (without Thomson and Rayleigh)
    #=
    α_cont = α_continuum.(line.λ0,
                          atmos.temperature,
                          atmos.electron_density*1.0,
                          LTE_pops[:, :, :, 1]*1.0,
                          LTE_pops[:, :, :, 3]*1.0)
    =#

    γ = γ_constant(line,
                   atmos.temperature,
                   (LTE_pops[:, :, :, 1] .+ LTE_pops[:, :, :, 2]),
                   atmos.electron_density)

    damping_λ = Array{Float64, 4}(undef, size(S_λ))
    for l in eachindex(line.λ)
        damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    profile = compute_voigt_profile(line, atmos, damping_λ, k)
    # profile = compute_doppler_profile(line, atmos, k)

    α_line = Array{Float64, 4}(undef, size(profile))u"m^-1"
    α_cont = copy(α_line)
    for l in eachindex(line.λ)
        α_line[l,:,:,:] = αline_λ(line,
                                  profile[l, :, :, :],
                                  populations[:, :, :, 2],
                                  populations[:, :, :, 1])

        α_cont[l,:,:,:] = α_absorption.(line.λ[l],
                                        atmos.temperature,
                                        atmos.electron_density*1.0,
                                        LTE_pops[:,:,:,1].+LTE_pops[:,:,:,2],
                                        LTE_pops[:,:,:,3]) .+
                          α_scattering.(line.λ[l],
                                        atmos.electron_density,
                                        LTE_pops[:,:,:,1])
    end
    α_tot = α_line .+ α_cont

    # plot_top_line(atmos, line, S_λ, α_tot, θ, ϕ, "LTE_line")
    for idλ in eachindex(line.λ)
        plot_top_intensity(atmos, line, S_λ, α_tot, θ, ϕ, idλ, "LTE_$idλ")
    end

end

DATA = "../data/bifrost_qs006023_s525.hdf5"
QUADRATURE = "../quadratures/ul7n12.dat"

compare(DATA, QUADRATURE);
# LTE_line(DATA)
print("")
