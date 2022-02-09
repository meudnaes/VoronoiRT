using Unitful
using HDF5

include("atmosphere.jl")

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

function write_to_file(populations::Array{<:NumberDensity, 4},
                       output_path::String)
    h5open(output_path, "r+") do file
        file["populations"][:,:,:,:] = ustrip.(populations .|> u"m^-3")
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

"""
    create_output_file(output_path::String, nλ::Int64, atmosphere_size::Tuple)
Initialise all output variables for the full atom mode.
"""
function create_output_file(output_path::String, nλ::Int64, atmosphere_size::Tuple,
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
    end
end
