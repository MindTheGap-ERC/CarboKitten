# ~/~ begin <<docs/src/denudation/empirical.md#src/Denudation/EmpiricalDenudationMod.jl>>[init]
module EmpiricalDenudationMod

import ..Abstract: DenudationType, denudation, redistribution
using ...Boxes: Box
using Unitful
export slope_kernel
# ~/~ begin <<docs/src/denudation/empirical.md#empirical-denudation>>[init]
function empirical_denudation(precip::Float64, slope::Any)
    local a = 9.1363
    local b = -0.008519
    local c = 580.51
    local d = 9.0156
    local e = -0.1245
    local f = 4.91086
    (a ./ (1 .+ exp.(b .* (precip .* 1000 .- c)))) .* (d ./ (1 .+ exp.(e .* (slope .- f)))) .* u"m/Myr"
end
# ~/~ end
# ~/~ begin <<docs/src/denudation/empirical.md#empirical-denudation>>[1]
@kwdef struct EmpiricalDenudation <: DenudationType
    precip::typeof(1.0u"m")
end
# ~/~ end
# ~/~ begin <<docs/src/denudation/empirical.md#empirical-denudation>>[2]
function slope_kernel(w::Any, cellsize::Float64)
    dzdx = (-w[1, 1] - 2 * w[2, 1] - w[3, 1] + w[1, 3] + 2 * w[2, 3] + w[3, 3]) / (8 * cellsize)
    dzdy = (-w[1, 1] - 2 * w[1, 2] - w[1, 3] + w[3, 1] + 2 * w[3, 2] + w[1, 1]) / (8 * cellsize)

    if abs(w[2, 2]) <= min.(abs.(w)...)
        return 0.0
    else
        atan(sqrt(dzdx^2 + dzdy^2)) * (180 / π)
    end
end
# ~/~ end

function denudation(::Box, p::EmpiricalDenudation, water_depth, slope, facies, state)
    precip = p.precip ./ u"m"
    denudation_rate = zeros(typeof(1.0u"m/Myr"), length(facies), size(slope)...)

    for idx in CartesianIndices(state.ca)
        f = state.ca[idx]
        if f == 0
            continue
        end
        if water_depth[idx] <= 0
            denudation_rate[f,idx] = empirical_denudation.(precip, slope[idx])
        end
    end
    return denudation_rate
end

function redistribution(::Box, p::EmpiricalDenudation, denudation_mass, water_depth)
    return nothing
end

end
# ~/~ end