using Plots
using UnitfulRecipes

include("line.jl")
include("functions.jl")
include("atmosphere.jl")
include("broadening.jl")
include("populations.jl")
include("characteristics.jl")

"""
    function circle_shape(x, y, r)

Function for a circle with centre coordinates (x, y) and radius r.
"""
function circle_shape(x, y, r)
    θ = LinRange(0, 2π, 500)
    x .+ r*cos.(θ), y .+ r*sin.(θ)
end

function plot_searchlight(k::Vector{Float64},
                          x,
                          y,
                          I,
                          R0,
                          title::String)

    @assert norm(k) ≈ 1 "Wrong direction vector!"

    RES = 300
    x_r = 0.5 - sign(k[1])*k[2]/k[1]
    if x_r < 0
        x_r = 1 - (ceil(x_r) - x_r)
    elseif x_r > 1
        x_r = x_r - floor(x_r)
    end

    y_r = 0.5 - sign(k[1])*k[3]/k[1]
    if y_r < 0
        y_r = 1 - (ceil(y_r) - y_r)
    elseif y_r > 1
        y_r = y_r - floor(y_r)
    end

    heatmap(ustrip.(x), ustrip.(y), ustrip.(transpose(I)),
            dpi=RES, title=title, xaxis="x", yaxis="y",
            aspect_ratio= :equal, clim=(0.,1.))
    plot!(circle_shape(ustrip(x_r), ustrip(y_r), ustrip(R0)),
          aspect_ratio = :equal,
          linecolor=:red,
          lw=2,
          label="")
    savefig("../img/compare_searchlight/$title")

end

function plot_top_intensity(atmos::Atmosphere,
                            S_λ::Array{<:UnitsIntensity_λ, 4},
                            α_tot::Array{<:PerLength, 4},
                            populations::Array{<:NumberDensity, 4},
                            θ::Float64,
                            ϕ::Float64,
                            idλ::Int,
                            title::String)

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    I_top = short_characteristics_up(k, S_λ[idλ, :, :, :], α_tot[idλ, :, :, :],
                                     atmos, I_0=S_λ[idλ, 1, :, :])

    I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[end, 2:end-1, 2:end-1]))

    heatmap(ustrip(atmos.x[2:end-1]),
         ustrip(atmos.y[2:end-1]),
         transpose(I_top),
         xaxis="x",
         yaxis="y",
         dpi=300,
         rightmargin=10Plots.mm,
         title=title,
         aspect_ratio=:equal)

    savefig("../img/compare_line/"*title)
end

function plot_top_line(atmos::Atmosphere,
                       line::HydrogenicLine,
                       S_λ::Array{<:UnitsIntensity_λ, 4},
                       α_tot::Array{<:PerLength, 4},
                       populations::Array{<:NumberDensity, 4},
                       θ::Float64,
                       ϕ::Float64,
                       title::String)

    nλ, nz, nx, ny = size(S_λ)
    I_λ = zero(S_λ)

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    for iλ in 1:nλ
        I_λ[iλ, :, :, :] = short_characteristics_up(k, S_λ[iλ, :, :, :], α_tot[iλ, :, :, :],
                                                    atmos, I_0=S_λ[iλ, 1, :, :])
    end

    start = 1
    stop = size(S_λ)[1]

    for idx in 1:5:size(S_λ)[end]
        for idy in 1:5:size(S_λ)[end-1]
            loc = "_"*string(idx)*"_"*string(idy)
            I_plot = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_λ))[:, end, idx, idy]
            scatter(ustrip.(line.λ)[start:stop],
                 ustrip.(I_plot)[start:stop],
                 xaxis="λ [nm]",
                 yaxis="Iλ [kW m^-2 nm^-1]",
                 dpi=300,
                 rightmargin=10Plots.mm,
                 title=title*loc)

            savefig("../img/compare_line/lines/"*title*loc)
        end
    end
end

"""
    read_quantities(DATA::String)

read quantities from simulation from a hdf5 file
"""
function read_quantities(DATA::String; periodic=true)
    atmos = Atmosphere(get_atmos(DATA; periodic=periodic, skip=1)...)
    local S_λ, populations
    h5open(DATA, "r") do file
        S_λ = read(file, "source_function")[:, :, :, :]*u"kW*m^-2*nm^-1"
        populations = read(file, "populations")[:, :, :, :]*u"m^-3"
    end

    if periodic
        S_λ_periodic = Array{Float64, 4}(undef, size(S_λ) .+ (0,0,2,2))u"kW*m^-2*nm^-1"
        populations_periodic = Array{Float64, 4}(undef, size(populations) .+ (0,2,2,0))u"m^-3"

        # fill inner box
        S_λ_periodic[:,:,2:end-1, 2:end-1] = S_λ
        # x-direction wall
        S_λ_periodic[:,:,1,2:end-1] = S_λ[:,:,end,:]
        S_λ_periodic[:,:,end,2:end-1] = S_λ[:,:,1,:]
        # y-direction wall
        S_λ_periodic[:,:,2:end-1,end] = S_λ[:,:,:,1]
        S_λ_periodic[:,:,2:end-1,1] = S_λ[:,:,:,end]
        # fix corners
        S_λ_periodic[:,:,1,1] .= S_λ[:,:,end,end]
        S_λ_periodic[:,:,1,end] .= S_λ[:,:,end,1]
        S_λ_periodic[:,:,end,1] .= S_λ[:,:,1,end]
        S_λ_periodic[:,:,end,end] .= S_λ[:,:,1,1]

        # fill inner box
        populations_periodic[:,2:end-1,2:end-1,:] = populations
        # x-direction wall
        populations_periodic[:,1,2:end-1,:] = populations[:,end,:,:]
        populations_periodic[:,end,2:end-1,:] = populations[:,1,:,:]
        # y-direction wall
        populations_periodic[:,2:end-1,end,:] = populations[:,:,1,:]
        populations_periodic[:,2:end-1,1,:] = populations[:,:,end,:]
        # fix corners
        populations_periodic[:,1,1,:] .= populations[:,end,end,:]
        populations_periodic[:,1,end,:] .= populations[:,end,1,:]
        populations_periodic[:,end,1,:] .= populations[:,1,end,:]
        populations_periodic[:,end,end,:] .= populations[:,1,1,:]

        return atmos, S_λ_periodic, populations_periodic
    else
        return atmos, S_λ, populations
    end
end

"""
    plotter(atmos::Atmosphere,
            S_λ::Array{<:UnitsIntensity_λ, 4},
            populations::Array{<:NumberDensity, 4},
            θ::Float64,
            ϕ::Float64)

plot simulations results
"""
function plotter(atmos::Atmosphere,
                 S_λ::Array{<:UnitsIntensity_λ, 4},
                 populations::Array{<:NumberDensity, 4},
                 θ::Float64,
                 ϕ::Float64,
                 title::String)

    idx = 12
    idy = 12

    line = HydrogenicLine(test_atom(50, 20)..., atmos.temperature)

    γ = γ_constant(line,
                   atmos.temperature,
                   (populations[:, :, :, 1].+populations[:, :, :, 2]),
                   atmos.electron_density)

    damping_λ = Array{Float64, 4}(undef, size(S_λ))
    for l in eachindex(line.λ)
        damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    profile = compute_voigt_profile(line, atmos, damping_λ, k)

    α_line = Array{Float64, 4}(undef, size(profile))u"m^-1"
    for l in eachindex(line.λ)
        α_line[l, :, :, :] = αline_λ(line,
                                     profile[l, :, :, :],
                                     populations[:, :, :, 1],
                                     populations[:, :, :, 2])
    end

    LTE_pops = LTE_populations(line, atmos)*1.0

    # Find continuum extinction and absorption extinction (only with Thomson and Rayleigh)
    α_cont = Array{Float64, 4}(undef, (length(line.λ), size(atmos.temperature)...))u"m^-1"
    for l in eachindex(line.λ)
        α_cont[l, :, :, :] = α_continuum.(line.λ[l],
                                          atmos.temperature*1.0,
                                          atmos.electron_density*1.0,
                                          LTE_pops[:, :, :, 1]*1.0,
                                          LTE_pops[:, :, :, 3]*1.0)
    end

    α_tot = α_line .+ α_cont

    plot_top_line(atmos, line, S_λ, α_tot, populations, θ, ϕ, title)

    # for i in 1:51
        # plot_top_intensity(atmos, S_λ, α_tot, populations, θ, ϕ, i, "i_map/top_intensity_regular"*string(i))
    # end
end

# plotter(read_quantities("../data/regular_line_1ray.h5", periodic=true)..., 10., 10., "Regular-Line")
plotter(read_quantities("../data/voronoi_line_1ray.h5", periodic=true)..., 10., 10., "Voronoi-Line")
