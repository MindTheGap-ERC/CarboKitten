# ~/~ begin <<docs/src/active-layer-transport.md#src/Components/ActiveLayer.jl>>[init]
@compose module ActiveLayer
@mixin WaterDepth, FaciesBase, SedimentBuffer

export disintegration, transportation

using ..Common
using CarboKitten.Transport.ActiveLayer: pde_stencil
using Unitful

@kwdef struct Facies <: AbstractFacies
    diffusion_coefficient::typeof(1.0u"m")
end

@kwdef struct Input <: AbstractInput
    disintegration_rate::Rate = 50.0u"m/Myr"
end

"""
    disintegration(input) -> f!

Prepares the disintegration step. Returns a function `f!(state::State)`. The returned function
modifies the state, popping sediment from the `sediment_buffer` and returns an array of `Amount`.
"""
function define_h(input::AbstractInput,state::AbstractState)
    max_h = input.disintegration_rate * input.time.Δt
    w = water_depth(input)(state)
    h = zeros(typeof(max_h), input.box.grid_size...)
    for i in CartesianIndices(input.box.grid_size)
        if w[i] > 0.0u"m"
            h[i] = min.(max_h, state.sediment_height[i])
        end
    end
    return h
end

function disintegration(input)
    #max_h = input.disintegration_rate * input.time.Δt
    output = Array{Float64, 3}(undef, n_facies(input), input.box.grid_size...)
        return function(state)
                #h = min.(max_h, state.sediment_height)
                h = define_h(input, state)
                state.sediment_height .-= h
                pop_sediment!(state.sediment_buffer, h ./ input.depositional_resolution .|> NoUnits, output)
                return output .* input.depositional_resolution 
        end
end

"""
    transportation(input::Input) -> f

Prepares the transportation step. Returns a function `f(state::State, active_layer)`,
transporting the active layer, returning a transported `Amount` of sediment.
"""
function transportation(input)
    x, y = box_axes(input.box)
    μ0 = input.initial_topography.(x, y')
    # We always return this array
    transported_output = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)
    stencils = [
        let stc = pde_stencil(input.box, f.diffusion_coefficient)
            (μ, p) -> @views stc(tuple.(μ, p[i,:,:]), transported_output[i,:,:])
        end for (i, f) in enumerate(input.facies) ]

    return function(state, active_layer::Array{Amount, 3})
        μ = state.sediment_height .+ μ0
        for stc in stencils
            stc(μ, active_layer)
        end

        return transported_output
    end
end

end
# ~/~ end
