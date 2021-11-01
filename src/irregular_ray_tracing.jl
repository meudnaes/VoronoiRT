include("voronoi_utils.jl")

function irregular_SC_up(sites::VoronoiSites, cells::AbstractArray{VoronoiCell},
                         tree::KDTree, I_0)
    # Traces rays through an irregular grid

    θ = 10
    ϕ = 170

    # start at the bottom
    # shoot rays through every grid cell

    # 1st layer
    z_upwind = sites.z_min

    # precalculate trigonometric functions
    cosθ = cos(θ)
    sinθ = sin(θ)

    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)

    for cell in cells
        # coordinate
        z = cell.z
        x = cell.x
        y = cell.y

        Δz = z - z_upwind

        # donwind position
        r = abs(Δz/cosθ)
        x_upwind = x + r*cosϕ*sinθ
        y_upwind = y + r*sinϕ*sinθ

        # number of natural neighbours
        k = cell.n

        # use the neighbours or the k closest points to interpolate
        neighbours = cell.neighbours


    end

end
