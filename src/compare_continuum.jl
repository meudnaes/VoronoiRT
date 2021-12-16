include("lambda_iteration.jl")

function compare(DATA, quadrature)

    function regular()
        global atmos
        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=6)...)

        maxiter = 10
        ϵ = 1e-5
        J_mean, S_λ, α_tot = Λ_regular(ϵ, maxiter, atmos, quadrature)

        return atmos
    end
    atmos = regular()

    function voronoi(atmos::Atmosphere)

        nx = length(atmos.x)
        nz = length(atmos.z)
        ny = length(atmos.y)

        n_sites = nz*nx*ny
        println(n_sites)
        positions = rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))

        sites_file = "../data/sites.txt"
        neighbours_file = "../data/neighbours.txt"
        println(positions)
        # write sites to file
        write_arrays(ustrip(positions[1, :]),
                     ustrip(positions[2, :]),
                     ustrip(positions[3, :]),
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
                             z_min, z_max,
                             x_min, x_max,
                             y_min, y_max,
                             n_sites)

        maxiter = 10
        ϵ = 1e-4
        J_mean, S_λ, α_tot = Λ_voronoi(ϵ, maxiter, atmos, quadrature)

    end
    voronoi(atmos)

end

compare("../data/bifrost_qs006023_s525_quarter.hdf5", "../quadratures/ul2n3.dat");
