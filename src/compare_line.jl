using Plots

include("io.jl")
include("line.jl")
include("functions.jl")
include("atmosphere.jl")
include("lambda_iteration.jl")

global my_seed = 2022
Random.seed!(my_seed)

function compare(DATA, quadrature)
    maxiter = 1
    ϵ = 1e-3

    θ = 10
    ϕ = 10

    n_skip = 1

    nλ_bb = 50
    nλ_bf = 20

    function regular()

        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=n_skip)...)

        global line
        line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., atmos.temperature)

        return

        REGULAR_DATA = "../data/regular_ul2n3_zero_radiation.h5"

        create_output_file(REGULAR_DATA, length(line.λ), size(atmos.temperature[:, 2:end-1, 2:end-1]), maxiter)
        write_to_file(atmos, REGULAR_DATA, ghost_cells=true)
        write_to_file(nλ_bb, "n_bb", REGULAR_DATA)
        write_to_file(nλ_bf, "n_bf", REGULAR_DATA)

        @time J_mean, S_λ, α_cont, populations = Λ_regular(ϵ, maxiter, atmos, line, quadrature, REGULAR_DATA)

        γ = γ_constant(line,
                       atmos.temperature,
                       (populations[:, :, :, 1] .+ populations[:, :, :, 2]),
                       atmos.electron_density)

        damping_λ = Array{Float64, 4}(undef, size(S_λ))
        for l in eachindex(line.λ)
            damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
        end

        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        profile = compute_voigt_profile(line, atmos, damping_λ, k)

        α_tot = Array{Float64, 4}(undef, size(profile))u"m^-1"
        for l in eachindex(line.λ)
            α_tot[l,:,:,:] = αline_λ(line,
                                        profile[l, :, :, :],
                                        populations[:, :, :, 1],
                                        populations[:, :, :, 2])
            α_tot[l,:,:,:] += α_cont
        end

        I_top = short_characteristics_up(k, S_λ[6,:,:,:], α_tot[6,:,:,:],
                                         atmos, I_0=S_λ[6,1,:,:])

        # plot_top_intensity(I_top, atmos.x, atmos.y, "regular_top")
        return 0
    end


    function voronoi()

        atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=n_skip)...)

        nx = length(atmos.x)
        nz = length(atmos.z)
        ny = length(atmos.y)

        n_sites = floor(Int, nz*nx*ny)

        positions = rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))

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
        # run(`./voro.sh $sites_file $neighbours_file
                       # $(x_min) $(x_max)
                       # $(y_min) $(y_max)
                       # $(z_min) $(z_max)`)

        # Voronoi grid
        sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions)...,
                             _initialise(positions, atmos)...,
                             z_min*1u"m", z_max*1u"m",
                             x_min*1u"m", x_max*1u"m",
                             y_min*1u"m", y_max*1u"m",
                             n_sites)

        line = HydrogenicLine(test_atom(nλ_bb, nλ_bf)..., sites.temperature)

        VORONOI_DATA = "../data/voronoi_ul2n3.h5"

        create_output_file(VORONOI_DATA, length(line.λ), n_sites, maxiter)
        write_to_file(nλ_bb, "n_bb", VORONOI_DATA)
        write_to_file(nλ_bf, "n_bf", VORONOI_DATA)
        write_to_file(sites, VORONOI_DATA)

        J_mean, S_λ, α_cont, populations = Λ_voronoi(ϵ, maxiter, sites, line, quadrature, VORONOI_DATA)

        γ = γ_constant(line,
                       sites.temperature,
                       (populations[:, 1].+populations[:, 2]),
                       sites.electron_density)

        damping_λ = Matrix{Float64}(undef, size(S_λ))
        for l in eachindex(line.λ)
            damping_λ[l, :] = damping.(γ, line.λ[l], line.ΔD)
        end

        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]

        profile = compute_voigt_profile(line, sites, damping_λ, k)

        α_tot = Matrix{Float64}(undef, size(profile))u"m^-1"
        for l in eachindex(line.λ)
            α_tot[l, :] = αline_λ(line,
                                  profile[l, :],
                                  populations[:, 1],
                                  populations[:, 2])
            α_tot[l,:] += α_cont
        end

        r_factor = 2
        atmos_from_voronoi, S_λ_grid, α_grid, populations_grid = Voronoi_to_Raster(sites, atmos, S_λ, α_tot, populations, r_factor)

        # plot_top_intensity(I_top, atmos_voronoi.x, atmos_voronoi.y, "irregular_top")
        return 0
    end

    # regular();
    voronoi();
end

DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
QUADRATURE = "../quadratures/ul2n3.dat"

compare(DATA, QUADRATURE);
print("")
