using Plots

include("io.jl")
include("line.jl")
include("functions.jl")
include("atmosphere.jl")
include("lambda_iteration.jl")

global my_seed = 2022
Random.seed!(my_seed)

function sample_from_extinction(atmos::Atmosphere,
                                λ0::Unitful.Length,
                                n_sites::Int)

    # Find continuum extinction and absorption extinction (without Thomson and Rayleigh)
    α_cont = α_continuum.(λ0,
                          atmos.temperature*1.0,
                          atmos.electron_density*1.0,
                          atmos.hydrogen_populations*1.0,
                          atmos.hydrogen_populations*1.0)

    positions = rejection_sampling(n_sites, atmos, ustrip.(α_cont))
    return positions
end

function compare(DATA, quadrature)
    maxiter = 100
    ϵ = 1e-3

    θ = 10
    ϕ = 10

    n_skip = 4

    nλ_bb = 50
    nλ_bf = 20

    function regular()

        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=n_skip)...)

        line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., atmos.temperature)

        REGULAR_DATA = "../data/regular_line_test-threads.h5"

        create_output_file(REGULAR_DATA, length(line.λ), size(atmos.temperature[:, 2:end-1, 2:end-1]), maxiter)
        write_to_file(nλ_bb, "n_bb", REGULAR_DATA)
        write_to_file(nλ_bf, "n_bf", REGULAR_DATA)

        J_mean, S_λ, α_cont, populations = Λ_regular(ϵ, maxiter, atmos, line, quadrature, REGULAR_DATA)

        γ = γ_constant(line,
                       atmos.temperature,
                       (populations[:, :, :, 1].+populations[:, :, :, 2]),
                       atmos.electron_density)

        damping_λ = Array{Float64, 4}(undef, size(S_λ))
        for l in eachindex(line.λ)
            damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
        end

        profile = compute_voigt_profile(line, atmos, damping_λ, θ*π/180, ϕ*π/180)

        α_tot = Array{Float64, 4}(undef, size(profile))u"m^-1"
        for l in eachindex(line.λ)
            α_tot[l,:,:,:] = αline_λ(line,
                                        profile[l, :, :, :],
                                        populations[:, :, :, 1],
                                        populations[:, :, :, 2])
            α_tot[l,:,:,:] += α_cont
        end

        I_top = short_characteristics_up(θ, ϕ, S_λ[6,:,:,:], α_tot[6,:,:,:],
                                         atmos, degrees=true, I_0=S_λ[6,1,:,:])

        # plot_top_intensity(I_top, atmos.x, atmos.y, "regular_top")



        write_to_file(populations[:, 2:end-1, 2:end-1, :], REGULAR_DATA)
        write_to_file(S_λ[:, :, 2:end-1, 2:end-1], REGULAR_DATA)
        write_to_file(atmos, REGULAR_DATA, ghost_cells=true)

        return 0
    end


    function voronoi()

        atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=n_skip)...)

        nx = length(atmos.x)
        nz = length(atmos.z)
        ny = length(atmos.y)

        n_sites = floor(Int, nz*nx*ny)

        global positions
        positions = sample_from_extinction(atmos, 121.0u"nm", n_sites)

        sites_file = "../data/sites_compare.txt"
        neighbours_file = "../data/neighbours_compare.txt"
        # write sites to file
        write_arrays(ustrip(positions[2, :]),
                     ustrip(positions[3, :]),
                     ustrip(positions[1, :]),
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

        line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., sites.temperature)

        VORONOI_DATA = "../data/voronoi_line_12.h5"
        create_output_file(VORONOI_DATA, length(line.λ), size(atmos_from_voronoi.temperature), maxiter)
        write_to_file(nλ_bb, "n_bb", VORONOI_DATA)
        write_to_file(nλ_bf, "n_bf", VORONOI_DATA)

        J_mean, S_λ, α_cont, populations = Λ_voronoi(ϵ, maxiter, sites, line, quadrature, VORONOI_DATA)

        γ = γ_constant(line,
                       sites.temperature,
                       (populations[:, 1].+populations[:, 2]),
                       sites.electron_density)

        damping_λ = Matrix{Float64}(undef, size(S_λ))
        for l in eachindex(line.λ)
            damping_λ[l, :] = damping.(γ, line.λ[l], line.ΔD)
        end

        profile = compute_voigt_profile(line, sites, damping_λ, θ*π/180, ϕ*π/180)

        α_tot = Matrix{Float64}(undef, size(profile))u"m^-1"
        for l in eachindex(line.λ)
            α_tot[l, :] = αline_λ(line,
                                  profile[l, :],
                                  populations[:, 1],
                                  populations[:, 2])
            α_tot[l,:] += α_cont
        end

        atmos_from_voronoi, S_λ_grid, α_grid, populations_grid = Voronoi_to_Raster(sites, atmos, S_λ, α_tot, populations, 3)

        # plot_top_intensity(I_top, atmos_voronoi.x, atmos_voronoi.y, "irregular_top")



        write_to_file(populations_grid, VORONOI_DATA)
        write_to_file(S_λ_grid, VORONOI_DATA)
        write_to_file(atmos_from_voronoi, VORONOI_DATA)

        return 0
    end

    regular();
    # voronoi();
end

DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
QUADRATURE = "../quadratures/ul7n12.dat"

compare(DATA, QUADRATURE);
print("")
