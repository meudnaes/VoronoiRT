"""
    function write_arrays(x::AbstractArray, y::AbstractArray, z::AbstractArray,
                          fname::String)

Writes the arrays z, x, and y to a file with filename fname.
Arrays are written in columns [ row number ] [ x ] [ y ] [ z ]
"""
function write_arrays(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64},
                      fname::String)

    if length(z) != length(y) || length(y) != length(x)
        println("Wrong shapes of input data")
        exit()
    end

    open(fname, "w") do io
        for i in 1:length(z)
            println(io, "$i\t$(x[i])\t$(y[i])\t$(z[i])")
        end
    end
end

function write_arrays(x::Vector{<:Unitful.Length}, y::Vector{<:Unitful.Length}, z::Vector{<:Unitful.Length},
                      fname::String)

    x = ustrip.(x)
    y = ustrip.(y)
    z = ustrip.(z)

    if length(z) != length(y) || length(y) != length(x)
        println("Wrong shapes of input data")
        exit()
    end

    open(fname, "w") do io
        for i in 1:length(z)
            println(io, "$i\t$(x[i])\t$(y[i])\t$(z[i])")
        end
    end
end

function write_boundaries(z_min, z_max, x_min, x_max, y_min, y_max, fname::String)
    open(fname, "w") do io
        println(io, "z_min = $z_min")
        println(io, "z_max = $z_max")
        println(io, "x_min = $x_min")
        println(io, "x_max = $x_max")
        println(io, "y_min = $y_min")
        println(io, "y_max = $y_max")
    end
end

"""
    write_to_file(Array,
                  output_path::String)
Write 'Array' to the output file.
"""
function write_to_file(S_λ::Array{<:UnitsIntensity_λ, 4},
                       output_path::String)
    h5open(output_path, "r+") do file
        file["source_function"][:,:,:,:] = ustrip.(S_λ .|> u"kW*m^-2*nm^-1")
    end
end
function write_to_file(S_λ::Matrix{<:UnitsIntensity_λ},
                       output_path::String)
    h5open(output_path, "r+") do file
        file["source_function"][:,:] = ustrip.(S_λ .|> u"kW*m^-2*nm^-1")
    end
end

function write_to_file(populations::Array{<:NumberDensity, 4},
                       output_path::String)
    h5open(output_path, "r+") do file
        file["populations"][:,:,:,:] = ustrip.(populations .|> u"m^-3")
    end
end
function write_to_file(populations::Matrix{<:NumberDensity},
                       output_path::String)
    h5open(output_path, "r+") do file
        file["populations"][:,:] = ustrip.(populations .|> u"m^-3")
    end
end

function write_to_file(atmos::Atmosphere,
                       output_path::String;
                       ghost_cells=false)
    h5open(output_path, "r+") do file
        for field in fieldnames(Atmosphere)
            if String(field) in ["z", "x", "y"]
                if ghost_cells && String(field) in ["x", "y"]
                    file[String(field)][:] = ustrip.(getfield(atmos, field))[2:end-1]
                else
                    file[String(field)][:] = ustrip.(getfield(atmos, field))
                end
            else
                if ghost_cells
                    file[String(field)][:, :, :] = ustrip.(getfield(atmos, field))[:, 2:end-1, 2:end-1]
                else
                    file[String(field)][:, :, :] = ustrip.(getfield(atmos, field))
                end
            end
        end
    end
end

function write_to_file(sites::VoronoiSites,
                       output_path::String)

    h5open(output_path, "r+") do file
        fieldnames = ["positions", "temperature", "electron_density",
                      "hydrogen_populations", "velocity_z", "velocity_x", "velocity_y"]

        for field in fieldnames
            if field == "positions"
                file[field][:, :] = ustrip.(getfield(sites, Symbol(field)))
            else
                file[field][:] = ustrip.(getfield(sites, Symbol(field)))
            end
        end


        file["boundaries"][:] = ustrip.([sites.z_min, sites.z_max,
                                         sites.x_min, sites.x_max,
                                         sites.y_min, sites.y_max])

    end

end


function write_to_file(difference::Float64,
                       iteration::Int,
                       output_path::String)
    h5open(output_path, "r+") do file
       file["convergence"][iteration] = difference
    end
end

function write_to_file(n::Int,
                       field::String,
                       output_path::String)
    h5open(output_path, "r+") do file
       file[field][:] = n
    end
end

function write_to_file(line::HydrogenicLine,
                       output_path::String)
    h5open(output_path, "r+") do file
        file["wavelength"][:] = ustrip.(line.λ .|> u"nm")
        file["line_center"][:] = ustrip(line.λ0 |> u"nm")
    end
end

"""
    create_output_file(output_path::String, nλ::Int64, atmosphere_size::Tuple)
Initialise all output variables for the full atom mode.
"""
function create_output_file(output_path::String, nλ::Int, atmosphere_size::Tuple,
                            maxiter::Int)

    nz, nx, ny = atmosphere_size

    h5open(output_path, "w") do file
        write(file, "source_function", Array{Float64, 4}(undef, (nλ, nz, nx, ny)))
        write(file, "populations", Array{Float64, 4}(undef, (nz, nx, ny, 3)))

        write(file, "z", Vector{Float64}(undef, nz))
        write(file, "x", Vector{Float64}(undef, nx))
        write(file, "y", Vector{Float64}(undef, ny))

        write(file, "temperature", Array{Float64, 3}(undef, (nz, nx, ny)))
        write(file, "hydrogen_populations", Array{Float64, 3}(undef, (nz, nx, ny)))
        write(file, "electron_density", Array{Float64, 3}(undef, (nz, nx, ny)))

        write(file, "velocity_z", Array{Float64, 3}(undef, (nz, nx, ny)))
        write(file, "velocity_x", Array{Float64, 3}(undef, (nz, nx, ny)))
        write(file, "velocity_y", Array{Float64, 3}(undef, (nz, nx, ny)))

        write(file, "convergence", zeros(Float64, maxiter+1))

        write(file, "n_bb", Vector{Int}(undef, 1))
        write(file, "n_bf", Vector{Int}(undef, 1))

        write(file, "wavelength", Vector{Float64}(undef, nλ))
        write(file, "line_center", Vector{Float64}(undef, 1))

        write(file, "time", Vector{Float64}(undef, 1))
    end
end

"""
    create_output_file(output_path::String, nλ::Int64, atmosphere_size::Tuple)
Initialise all output variables for the full atom mode.
"""
function create_output_file(output_path::String, nλ::Int, n_sites::Int,
                            maxiter::Int)

    h5open(output_path, "w") do file
        write(file, "source_function", Matrix{Float64}(undef, (nλ, n_sites)))
        write(file, "populations", Matrix{Float64}(undef, (n_sites, 3)))

        write(file, "positions", Matrix{Float64}(undef, (3, n_sites)))

        write(file, "temperature", Vector{Float64}(undef, n_sites))
        write(file, "hydrogen_populations", Vector{Float64}(undef, n_sites))
        write(file, "electron_density", Vector{Float64}(undef, n_sites))

        write(file, "velocity_z", Vector{Float64}(undef, n_sites))
        write(file, "velocity_x", Vector{Float64}(undef, n_sites))
        write(file, "velocity_y", Vector{Float64}(undef, n_sites))

        write(file, "boundaries", Vector{Float64}(undef, 6))

        write(file, "convergence", zeros(Float64, maxiter+1))

        write(file, "n_bb", Vector{Int}(undef, 1))
        write(file, "n_bf", Vector{Int}(undef, 1))

        write(file, "wavelength", Vector{Float64}(undef, nλ))
        write(file, "line_center", Vector{Float64}(undef, 1))

        write(file, "time", Vector{Float64}(undef, 1))
    end
end
