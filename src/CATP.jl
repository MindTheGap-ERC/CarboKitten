# ~/~ begin <<docs/src/submarine-transport.md#src/CATP.jl>>[init]
module CATP

using ..Vectors
using ..Stencil: Periodic
using ..Utility
using ..Stencil
using ..Burgess2013.CA
using ..SedimentStack
using ..Transport
using ..BoundaryTrait

using HDF5
using Printf

using .Iterators: drop, peel, partition, map, take, product
# ~/~ begin <<docs/src/submarine-transport.md#cat-facies>>[init]
struct Facies
  viability_range::Tuple{Int,Int}
  activation_range::Tuple{Int,Int}

  maximum_growth_rate::Float64
  extinction_coefficient::Float64
  saturation_intensity::Float64

  grain_size::Float64
  excess_density::Float64
  # critical_angle::Float64

  # magnitude of critical stress, should be ΔρDg sin α where α is the critical angle
  critical_stress::Float64
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-input>>[init]
@kwdef struct Input
  sea_level
  subsidence_rate
  initial_depth

  grid_size::NTuple{2,Int}
  boundary::Type
  phys_scale::Float64     # km / pixel
  Δt::Float64             # Ma / timestep ?

  time_steps::Int
  write_interval::Int

  facies::Vector{Facies}
  insolation::Float64

  Δz::Float64
  buffer_depth::Int
  # Function of depth
  disintegration_rate
  # wave shear stress function, is a function of depth, and should return a Vec2
  wave_shear_stress
  # Gravitational accelleration
  g::Float64
  # for the transport step, parcels of sediment are converted into particles
  # by setting this parameter to something >1 you can control how many particles
  # are created on each axis. So number of particles is grid width * height * transport_subsample^2.
  transport_subsample::Int
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-state>>[init]
mutable struct State
  time::Float64
  height::Array{Float64,2}            # x y
  sediment::Array{Float64,4}          # x y z (f..., t)
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-frame>>[init]
struct ProductFrame
  production::Array{Float64,3}        # facies x y
end

Base.:+(a::ProductFrame, b::ProductFrame) = ProductFrame(a.production .+ b.production)

# struct Properties
#   excess_density::Float64
#   grain_size::Float64
# end

const Particle = Transport.Particle{Nothing}
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-frame>>[1]
struct ReductionFrame
  disintegration::Array{Float64,2}   # x y
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-update>>[init]
function deposit_material(
    grid_size::NTuple{2, Int},
    Δz::Float64,
    s::State,
    material::AbstractArray{Float64,3})

  Threads.@threads for idx in CartesianIndices(grid_size)
    production = material[Tuple(idx)..., :] .* (1.0 / Δz)
    push_sediment!(s.sediment, production)
    s.height[idx] .-= sum(material[Tuple(idx)..., :]) * Δt
  end
end

function remove_material(
    grid_size::NTuple{2, Int},
    n_facies::Int,
    Δz::Float64,
    s::State,
    thickness::AbstractArray{Float64,2})

  material = Array{Float64,3}(undef, grid_size..., n_facies)
  Threads.@threads for idx in CartesianIndices(grid_size)
    material[Tuple(idx)..., :] = pop_sediment!(s.sediment, thickness[idx])
  end
  return material
end

function deposit_updater(input::Input)
  function (s::State, Δ::ProductFrame)
    deposit_material(input.grid_size, input.Δz, s, Δ.production)
  end
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-update>>[1]
function remove_updater(input::Input)
  function (s::State, Δ::ReductionFrame)
    loose = remove_material(input.grid_size, length(input.facies), input.Δz, s, Δ.disintegration)
    # need physics to disperse loose material
    return ProductFrame(loose)
  end
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-update>>[2]
function time_updater(input::Input)
  function (s::State, Δ::Nothing = nothing)
    s.height .+= input.subsidence_rate * input.Δt
    s.time += input.Δt
  end
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-propagator>>[init]
function particles(input::Input, Δ::ProductFrame)
  Channel{Particle}() do ch
    for (i, (idx, mass)) in enumerate(pairs(Δ.production))
      subgrid_spacing = 1.0 / input.transport_subsample
      subgrid_axis = 0.0:subgrid_spacing:1.0 - subgrid_spacing
      subgrid = product(subgrid_axis, subgrid_axis)
      for (j, dx) in enumerate(subgrid)
        facies_type = idx[1]
        p = (x=(idx[2]-1 + dx[1]) * input.phys_scale, y=(idx[3]-1 + dx[2]) * input.phys_scale)
        θ = input.facies[facies_type].critical_stress
        put!(ch, Particle(p, mass * subgrid_spacing^2, θ, facies_type, nothing))
      end
    end
  end
end

function stress(input::Input, s::State)
  phys_size = (x = input.grid_size[1] * input.phys_scale,
               y = input.grid_size[2] * input.phys_scale)
  box = Transport.Box(input.grid_size, phys_size, input.phys_scale)

  function (p::Particle)
    z, ∇ = Transport.interpolate(input.boundary, box, s.height, p.position)
    α = atan(abs(∇))
    Ĝ = -∇ / abs(∇)

    τ_wave = input.wave_shear_stress(z)

    Δρ = input.facies[p.facies].excess_density
    D = input.facies[p.facies].grain_size
    g = input.g
    τ_grav = (Δρ * D * g * sin(α)) * Ĝ

    return τ_grav + τ_wave
  end
end

function submarine_transport(input::Input)
  # This function does not modify the state, rather it transports the
  # sediments in a given product frame and gives a new product frame
  # as output.
  box = Transport.Box(input.grid_size, (
    x=input.grid_size[1] * input.phys_scale,
    y=input.grid_size[2] * input.phys_scale), input.phys_scale)
  function (s::State, Δ::ProductFrame)  # -> ProductFrame
    output = zeros(Float64, size(Δ.production)...)
    transport = Transport.transport(input.boundary, box, stress(input, s))
    deposit = Transport.deposit(input.boundary, box, output)

    for p in particles(input, Δ)
      p |> transport |> deposit
    end

    return ProductFrame(output)
  end
end
# ~/~ end
function production_propagator(input::Input)
    # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-init-propagator>>[init]
    n_facies = length(input.facies)
    ca_init = rand(0:n_facies, input.grid_size...)
    ca = drop(run_ca(Periodic{2}, input.facies, ca_init, 3), 20)

    function water_depth(s::State)
        s.height .- input.sea_level(s.time)
    end
    # ~/~ end
    function (s::State)  # -> Frame
        # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-propagate>>[init]
        result = zeros(Float64, input.grid_size..., n_facies)
        facies_map, ca = peel(ca)
        w = water_depth(s)
        Threads.@threads for idx in CartesianIndices(facies_map)
            f = facies_map[idx]
            if f == 0
                continue
            end
            result[Tuple(idx)..., f] = production_rate(input.insolation, input.facies[f], w[idx])
        end
        return Frame(result)
        # ~/~ end
    end
end

function initial_state(input::Input)  # -> State
    height = zeros(Float64, input.grid_size...)
    for i in CartesianIndices(height)
        height[i] = input.initial_depth(i[2] * input.phys_scale)
    end
    n_facies = length(input.facies)
    return State(0.0, height, zeros(Float64, input.grid_size..., input.buffer_depth, n_facies))
end

function disintegration_propagator(input::Input)
  function (s::State)
    return ReductionFrame(input.disintegration_rate.(s.height))
  end
end

struct Snapshot
  state::State
  removed::Array{Float64,2}
  deposited::Array{Float64,3}
end

function run_model(input::Input)
  transport = submarine_transport(input)
  p_produce = production_propagator(input)
  p_disintegrate = disintegration_propagator(input)
  state = initial_state(input)
  u_time = time_updater(input)
  u_remove = remove_updater(input)
  u_deposit = deposit_updater(input)

  Channel{Snapshot}() do ch
    while True
      Δ_produced = p_produce(state)
      reduction = p_disintegrate(state)
      Δ_removed = u_remove(state, reduction)
      Δ_transported = transport(state, Δ_produced + Δ_removed)
      u_deposit(state, Δ_transported)
      u_time(state)
      put!(ch, Snapshot(state, Δ_removed, Δ_transported))
    end
  end
end

function main(input::Input, output::String)
    x_axis = (0:(input.grid_size[2]-1)) .* input.phys_scale
    y_axis = (0:(input.grid_size[1]-1)) .* input.phys_scale
    initial_height = input.initial_depth.(x_axis)
    n_writes = input.time_steps

    h5open(output, "w") do fid
        gid = create_group(fid, "input")
        gid["x"] = collect(x_axis)
        gid["y"] = collect(y_axis)
        gid["height"] = collect(initial_height)
        gid["t"] = collect((0:(n_writes-1)) .* (input.Δt * input.write_interval))
        attr = attributes(gid)
        attr["delta_t"] = input.Δt
        attr["time_steps"] = input.time_steps
        attr["subsidence_rate"] = input.subsidence_rate

        results = run_model(input)
        for (step, snapshot) in enumerate(take(results, n_writes))
          gid = create_group(fid, @sprintf "%010u" step)
          attr = attributes(gid)
          attr["time"] = snapshot.state.time
          gid["height"] = snapshot.state.height
          gid["buffer"] = snapshot.state.sediment
          gid["removed"] = snapshot.removed
          gid["deposited"] = snapshot.deposited
        end
    end
end

end  # CATP
# ~/~ end