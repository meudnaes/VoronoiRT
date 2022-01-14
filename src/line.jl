using Transparency
include("functions.jl")

struct Line
    u::Int
    l::Int
    lineData::AtomicLine
    doppler_width::Array{<:Unitful.Length, 3}        # (n_lines, nz, nx, ny)
    damping_constant::Array{<:PerArea,3}             # (n_lines, nx, ny)
end
