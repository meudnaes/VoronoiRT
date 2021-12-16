include("lambda_iteration.jl")

function compare(DATA, quadrature)

    function regular()
        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=2)...)

        maxiter = 1000
        ϵ = 1e-5
        J_mean, S_λ, α_tot = Λ_regular(ϵ, maxiter, atmos, quadrature)

        return atmos
    end
    atmos = regular()

    function voronoi(atmos::Atmosphere)

        maxiter = 1000
        ϵ = 1e-4
        J_mean, S_λ, α_tot = Λ_voronoi(ϵ, maxiter, atmos, quadrature)

    end
    #voronoi(atmos)

end

compare("../data/bifrost_qs006023_s525_quarter.hdf5", "../quadratures/ul2n3.dat");
