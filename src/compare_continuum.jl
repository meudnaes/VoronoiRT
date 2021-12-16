include("characteristics.jl")
include("irregular_ray_tracing.jl")

function compare()

    function regular()
        DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=4)...)
    end
    regular()

end

compare()
