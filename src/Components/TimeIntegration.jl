# ~/~ begin <<docs/src/components/time.md#src/Components/TimeIntegration.jl>>[init]
@compose module TimeIntegration
    using ..Common

    @kwdef struct Input <: AbstractInput
        time::TimeProperties
    end

    mutable struct State <: AbstractState
        time::Time
    end

    State(input::AbstractInput) = State(0.0u"Myr")
end
# ~/~ end