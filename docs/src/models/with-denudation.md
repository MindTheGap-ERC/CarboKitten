# Model including denudation

This is largely identical to the ALCAP model.

``` {.julia file=src/Models/WithDenudation.jl}
@compose module WithDenudation
@mixin Tag, H5Writer, CAProduction, ActiveLayer, Denudation

using ..Common
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth
using ModuleMixins: @for_each
using ...Stencil
using ...BoundaryTrait
using ...Denudation.EmpiricalDenudationMod: slope_kernel
export Input, Facies

function initial_state(input::Input)
    ca_state = CellularAutomaton.initial_state(input)
    for _ in 1:20
        CellularAutomaton.step!(input)(ca_state)
    end

    sediment_height = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)

    return State(
        step=0, sediment_height=sediment_height,
        sediment_buffer=sediment_buffer,
        ca=ca_state.ca, ca_priority=ca_state.ca_priority)
end

function step!(input::Input)
    step_ca! = CellularAutomaton.step!(input)
    disintegrate! = disintegration(input)
    produce = production(input)
    transport = transportation(input)
    denudate = denudation(input)
    redistribute = redistribution(input)
    water_depth_fn = water_depth(input)
    slopefn = slope_function(input,input.box)
    # Somehow deal with the units here


    #slopefn = stencil(Float64, BT, (3, 3), slope_kernel)
    slope = Array{Float64}(undef, input.box.grid_size...)
    denuded_sediment = Array{Float64}(undef, n_facies(input), input.box.grid_size...)

    function (state::State)
        if mod(state.step, input.ca_interval) == 0
            step_ca!(state)
        end

        w = water_depth_fn(state) ./u"m"
        slopefn(w, slope, input.box.phys_scale ./ u"m")

        # submarine: production and transport
        p = produce(state)
        d = disintegrate!(state)

        active_layer = p .+ d
        sediment = transport(state, active_layer)


        # subaerial: denudation and redistribution
        # this code should go into the Denudation component
        #min.(sum(denudate(state,w,slope),dims=1), state.sediment_height)
        denudation_mass = denudate(state,w,slope)
        if denudation_mass !== nothing
            denudation_mass = denudate(state,w,slope) |> x -> sum(x,dims=1) |> x -> dropdims(x,dims=1) |> x -> min.(x, state.sediment_height)

            state.sediment_height -= denudation_mass
            pop_sediment!(state.sediment_buffer, denudation_mass ./ input.depositional_resolution .|> NoUnits, denuded_sediment)

            d += denuded_sediment .* input.depositional_resolution

            redistribution_mass = redistribute(state, w, denuded_sediment .* input.depositional_resolution)
            if redistribution_mass !== nothing
                sediment .+= redistribution_mass
            end
        end

        push_sediment!(state.sediment_buffer, sediment ./ input.depositional_resolution .|> NoUnits)

        state.sediment_height .+= sum(sediment; dims=1)[1,:,:]
        state.step += 1

        return H5Writer.DataFrame(
            production = p,
            disintegration = d,
            deposition = sediment)
    end
end

function write_header(fid, input::AbstractInput)
    @for_each(P -> P.write_header(fid, input), PARENTS)
end
end
```
