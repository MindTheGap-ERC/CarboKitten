module Abstract

using ...BoundaryTrait: Boundary

abstract type DenudationType end

"""
    denudation(param, state)

Computes the denudation for a single time-step, given denudation parameters
`param` and a simulation state `state`. `param` should have a `DenudationType`
type and `state` should contain the `height` property and `sealevel`.
"""
function denudation(::Type{BT}, param::DenudationType, water_depth, slope, facies) where {BT<:Boundary}
    error("Abstract `denudation` function called.")
end

function redistribution(::Type{BT}, param::DenudationType, water_depth, slope, facies) where {BT<:Boundary}
    error("Abstract `redistribution` function called.")
end

end  # module
