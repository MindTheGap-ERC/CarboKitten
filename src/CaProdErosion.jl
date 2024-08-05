# ~/~ begin <<docs/src/ca-prod-with-erosion.md#src/CaProdErosion.jl>>[init]
module CaProdErosion

using CarboKitten
using ..Stencil: Periodic, stencil
using ..Utility
using ..BoundaryTrait: Boundary
#using ..BS92: sealevel_curve
using ..Denudation: denudation, calculate_redistribution
using ..Denudation.EmpericalDenudation: slope_kernel
using ..Burgess2013
using ..Burgess2013.CA: step_ca, run_ca
using ..InputConfig: Input, DenudationType
using ..Config: Box
using Unitful
using HDF5
using .Iterators: drop, peel, partition, map, take
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-input>>[init]
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-frame>>[init]

abstract type Frame end
struct ProductionFrame <: Frame
    production::Array{typeof(1.0u"m/Myr"),3}
end

struct DenudationFrame <: Frame
    denudation::Array{typeof(1.0u"m/Myr"),2}
    redistribution::Array{typeof(1.0u"m/Myr"),2}
end

struct OutputFrame 
    production::Array{typeof(1.0u"m/Myr"),3}
    denudation::Array{typeof(1.0u"m/Myr"),2}
    redistribution::Array{typeof(1.0u"m/Myr"),2}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-state>>[init]
# FIXME: deduplicate
mutable struct State
    time::typeof(1.0u"Myr")
    ca::Array{Int}
    ca_priority::Vector{Int}
    height::Array{typeof(1.0u"m"),2}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-model>>[init]
# FIXME: deduplicate
function initial_state(input::Input)  # -> State
    height = zeros(Float64, input.box.grid_size...) * u"m"
    for i in CartesianIndices(height)
        height[i] = input.initial_depth(i[2] * input.box.phys_scale) 
    end
    n_facies = length(input.facies)
    ca = rand(0:n_facies, input.box.grid_size...)
    state = State(0.0u"Myr", ca, 1:n_facies, height)

    step = step_ca(input.box, input.facies)
    for _ = 1:20
        step(state)
    end
    return state
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-model>>[1]

# propagator for production
function prod_propagator(input::Input,box::Box{BT}) where {BT<:Boundary}
    n_facies = length(input.facies)

    function water_depth(s::State)
        sea_level = input.sea_level(s.time) .* u"m"
        s.height .- sea_level
    end

    function(s::State)
    production = zeros(typeof(0.0u"m/Myr"), box.grid_size..., n_facies)
    w = water_depth(s)
    for idx in CartesianIndices(s.ca)
        f = s.ca[idx]
        if f == 0
            continue
        end
        production[Tuple(idx)..., f] = production_rate(input.insolation, input.facies[f], w[idx])
    end
    return ProductionFrame(production)
    end

    
end


function denu_propagator(input::Input, box::Box{BT}) where {BT <: Boundary}

    function water_depth(s::State)
        sea_level = input.sea_level(s.time) .* u"m"
        s.height .- sea_level
    end

    function get_inf_map(s::State,input::Input)
        w = water_depth(s) ./ u"m"
        inf_map = ones(size(w)...)
        for idx in CartesianIndices(s.ca)
            f = s.ca[idx]
            if f == 0
                continue
            end
            inf_map[idx] = input.facies[f].infiltration_coefficient
        end
        return inf_map
    end

    function (s::State)
        w = water_depth(s) ./ u"m"
        slope = zeros(Float64, box.grid_size...)
        slopefn = stencil(Float64, BT, (3, 3), slope_kernel)
        slopefn(w, slope, box.phys_scale ./u"m")
        denudation_mass = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
        w = water_depth(s) ./ u"m"

        for idx in CartesianIndices(s.ca)
            f = s.ca[idx]
            if f == 0
                continue
            end

            if w[idx] >= 0
            (denudation_mass[idx]) = denudation(box, input.denudationparam, w[idx], slope[idx],input.facies[f])
            end
        end
    
        inf_map = get_inf_map(s,input)
        redistribution_mass = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
        (redistribution_mass) = calculate_redistribution(box,input.denudationparam,w,slope,inf_map)
    
    return DenudationFrame(denudation_mass,redistribution_mass)
    
    end
    
end

# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-model>>[2]
function updater(input)
    n_facies = length(input.facies)
    function update(state, Δ::ProductionFrame)
        state.height .-= sum(Δ.production; dims=3) .* input.time.Δt
        state.time += input.time.Δt
    end

    function update(state, Δ::DenudationFrame)
        # FIXME: implement
        state.height .+= Δ.denudation .* input.time.Δt
        state.height .-= Δ.redistribution .* input.time.Δt
        state.time += input.time.Δt
    end

    update
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-model>>[3]
function run_model(input::Input,box::Box{BT}) where {BT <:Boundary}
    Channel{OutputFrame}() do ch
        s = initial_state(input)
        p = prod_propagator(input,box)
        d = denu_propagator(input,box)  # FIXME: implement
        u = updater(input)

        while true
            Δ_prod = p(s)
            u(s, Δ_prod)
            Δ_denu = d(s)
            u(s, Δ_denu)
            put!(ch, OutputFrame(Δ_prod.production, Δ_denu.denudation, Δ_denu.redistribution))
        end
    end
end
# ~/~ end

function stack_frames(fs::Vector{OutputFrame})  # -> Frame
    OutputFrame(sum(f.production for f in fs),sum(f.denudation for f in fs),sum(f.redistribution for f in fs))#
end

function main(input::Input, output::String)
    x_axis = (0:(input.box.grid_size[2]-1)) .* input.box.phys_scale
    y_axis = (0:(input.box.grid_size[1]-1)) .* input.box.phys_scale
    initial_height = input.initial_depth.(x_axis)
    n_writes = input.time.steps ÷ input.time.write_interval
    t = collect((0:(n_writes-1)) .* (input.time.Δt * input.time.write_interval))

    h5open(output, "w") do fid
        gid = create_group(fid, "input")
        gid["x"] = collect(x_axis) |> in_units_of(u"m")
        gid["y"] = collect(y_axis) |> in_units_of(u"m")
        gid["height"] = collect(initial_height) |> in_units_of(u"m")
        gid["t"] =  t |> in_units_of(u"Myr")
        attr = attributes(gid)
        attr["delta_t"] = input.time.Δt |> in_units_of(u"Myr")
        attr["write_interval"] = input.time.write_interval
        attr["time_steps"] = input.time.steps
        attr["sea_level"] = input.sea_level.(t) 
        attr["subsidence_rate"] = input.subsidence_rate |> in_units_of(u"m/Myr")
        println("Subsidence rate saved successfully.")


        n_facies = length(input.facies)
        ds = create_dataset(fid, "sediment", datatype(Float64),
            dataspace(input.box.grid_size..., n_facies, input.time.steps),
            chunk=(input.box.grid_size..., n_facies, 1))
        denudation = create_dataset(fid, "denudation", datatype(Float64),
            dataspace(input.box.grid_size..., input.time.steps),
           chunk=(input.box.grid_size..., 1))
        redistribution = create_dataset(fid, "redistribution", datatype(Float64),
           dataspace(input.box.grid_size..., input.time.steps),
          chunk=(input.box.grid_size..., 1))
        box = input.box

        results = map(stack_frames, partition(run_model(input,box), input.time.write_interval))
        for (step, frame) in enumerate(take(results, n_writes))
            ds[:, :, :, step] = frame.production |> in_units_of(u"m/Myr")  
            denudation[:,:,step] = frame.denudation |> in_units_of(u"m/kyr")
            redistribution[:,:,step] = frame.redistribution |> in_units_of(u"m/kyr")
        end
    end
end

end # CaProd
# ~/~ end