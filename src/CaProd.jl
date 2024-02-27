# ~/~ begin <<docs/src/ca-with-production.md#src/CaProd.jl>>[init]
module CaProd

using CarboKitten
using ..Stencil: Periodic, stencil
using ..Utility
#using ..BS92: sealevel_curve
using ..Stencil
using ..Burgess2013
using ..EmpericalDenudation
using ..CarbDissolution
using PhysicalErosion

using HDF5
using .Iterators: drop, peel, partition, map, take

# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-input>>[init]
@kwdef struct Input
    sea_level
    subsidence_rate
    initial_depth

    grid_size::NTuple{2,Int}
    phys_scale::Float64
    Δt::Float64

    time_steps::Int
    write_interval::Int

    facies::Vector{Facies}
    insolation::Float64

    temp #temperature
    precip #precipitation
    pco2 #co2
    alpha #reaction rate

end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-frame>>[init]
struct Frame
    production::Array{Float64,3}
    denudation::Array{Float64,2}
    redistribution::Array{Float64,2}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-state>>[init]
mutable struct State
    time::Float64
    height::Array{Float64,2}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[init]
function initial_state(input::Input)  # -> State
    height = zeros(Float64, input.grid_size...)
    for i in CartesianIndices(height)
        height[i] = input.initial_depth(i[2] * input.phys_scale)
    end
    return State(0.0, height)
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[1]
function propagator(input::Input)
    # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-init-propagator>>[init]
    n_facies = length(input.facies)
    ca_init = rand(0:n_facies, input.grid_size...)
    ca = drop(run_ca(Periodic{2}, input.facies, ca_init, 3), 20)

    function water_depth(s::State)
        s.height .- input.sea_level(s.time)
    end
    # prepare functions for erosion
    # ~/~ end
    slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel)
    function (s::State)  # -> Frame
        # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-propagate>>[init]
        production = zeros(Float64, input.grid_size..., n_facies)
        denudation = zeros(Float64, input.grid_size...)
        redistribution = zeros(Float64, input.grid_size...)
        slope = zeros(Float64, input.grid_size...)
        facies_map, ca = peel(ca)
        w = water_depth(s)
        slopefn(w,slope,input.phys_scale) # slope is calculated with square so no need for -w
        redis = mass_erosion(Float64,Periodic{2},slope,(3,3))
        Threads.@threads for idx in CartesianIndices(facies_map)
            f = facies_map[idx]
                if f == 0
                    continue
                end
            if w[idx] > 0.0
                production[Tuple(idx)..., f] = production_rate(input.insolation, input.facies[f], w[idx])
            else
                #denudation[Tuple(idx)...] = dissolution(input.temp,input.precip,input.alpha,input.pco2,w[idx],input.facies[f])
                denudation[Tuple(idx)...] = emperical_denudation(input.precip, slope[idx])
                denudation[Tuple(idx)...] = physical_erosion(slope[idx],input.facies.inf)
                redistribution[Tuple(idx)...] = total_mass_redistribution(redis, slope)
            end
        end
        return Frame(production, denudation, redistribution)#
        # ~/~ end
    end
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[2]
function updater(input::Input)
    n_facies = length(input.facies)
    function (s::State, Δ::Frame)
        s.height .-= sum(Δ.production; dims=3) .* input.Δt
        s.height .+= Δ.denudation .* input.Δt  #number already in kyr
        s.height ._= Δ.redistribution .* input.Δt
        s.height .+= input.subsidence_rate * input.Δt
        s.time += input.Δt
    end
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[3]
function run_model(input::Input)
    Channel{Frame}() do ch
        s = initial_state(input)
        p = propagator(input)
        u = updater(input)

        while true
            Δ = p(s)
            put!(ch, Δ)
            u(s, Δ)
        end
    end
end
# ~/~ end

function stack_frames(fs::Vector{Frame})  # -> Frame
    Frame(sum(f.production for f in fs),sum(f.denudation for f in fs),sum(f.wredistribution for f in fs))#
end

function main(input::Input, output::String)
    x_axis = (0:(input.grid_size[2]-1)) .* input.phys_scale
    y_axis = (0:(input.grid_size[1]-1)) .* input.phys_scale
    initial_height = input.initial_depth.(x_axis)
    n_writes = input.time_steps ÷ input.write_interval

    h5open(output, "w") do fid
        gid = create_group(fid, "input")
        gid["x"] = collect(x_axis)
        gid["y"] = collect(y_axis)
        gid["height"] = collect(initial_height)
        gid["t"] = collect((0:(n_writes-1)) .* (input.Δt * input.write_interval))
        attr = attributes(gid)
        attr["delta_t"] = input.Δt
        attr["write_interval"] = input.write_interval
        attr["time_steps"] = input.time_steps
        attr["subsidence_rate"] = input.subsidence_rate

        n_facies = length(input.facies)
        ds = create_dataset(fid, "sediment", datatype(Float64),
            dataspace(input.grid_size..., n_facies, input.time_steps),
            chunk=(input.grid_size..., n_facies, 1))
        erosion = create_dataset(fid, "erosion", datatype(Float64),
            dataspace(input.grid_size..., input.time_steps),
           chunk=(input.grid_size..., 1))

        results = map(stack_frames, partition(run_model(input), input.write_interval))
        for (step, frame) in enumerate(take(results, n_writes))
            ds[:, :, :, step] = frame.production
            erosion[:,:,step] = frame.erosion
        end
    end
end

end # CaProd
# ~/~ end