# ~/~ begin <<docs/src/denudation/chemical.md#src/Denudation/DissolutionMod.jl>>[init]
# based on Kaufman 2002, Geomorphology
module DissolutionMod

import ..Abstract: DenudationType, denudation, redistribution
using ...BoundaryTrait: Boundary
using ...Config: Box

using Unitful

# TODO add units

@kwdef struct Dissolution <: DenudationType
    temp::typeof(1.0u"K")
    precip::typeof(1.0u"m")
    pco2::typeof(1.0u"atm")
    reactionrate::typeof(1.0u"m/yr")
end

# Kaufmann 2002, Table 2
function karst_denudation_parameters(temp::Float64)
    A = -0.4883 + 8.074 * 0.0001 * (temp - 273.0)
    B = -0.3241 + 1.6 * 0.0001 * (temp - 273.0)
    IA = 0.1 # ion activity

    (K1=10^(-356.3094 - 0.06091964 * temp + 21834.37 / temp + 126.8339 * log10(temp) - 1684915 / (temp^2)),
        K2=10^(-107.881 - 0.03252849 * temp + 5151.79 / temp + 38.92561 * log10(temp) - 563713.9 / (temp^2)),
        KC=10^(-171.9065 - 0.077993 * temp + 2839.319 / temp + 71.595 * log10(temp)),
        KH=10^(108.3865 + 0.01985076 * temp - 6919.53 / temp - 40.4515 * log10(temp) + 669365 / (temp^2)),
        activity_Ca=10^(-4A * sqrt(IA) / (1 + 10^(-8) * B * sqrt(IA))),
        activity_Alk=10^(-A * sqrt(IA) / (1 + 5.4 * 10^(-8) * B * sqrt(IA))))
end

#calculate ceq and Deq, Kaufman 2002
function equilibrium(temp::Float64, pco2::Float64, precip::Float64, facies)
    p = karst_denudation_parameters(temp)
    mass_density = facies.mass_density ./ u"kg/m^3"
    eq_c = (pco2 .* (p.K1 * p.KC * p.KH) ./ (4 * p.K2 * p.activity_Ca .* (p.activity_Alk)^2)) .^ (1 / 3) 
    eq_d = 1000 * precip .* facies.infiltration_coefficient * 40 * 1000 .* eq_c ./ mass_density 
    (concentration=eq_c, denudation=eq_d)
end

function dissolution(temp, precip, pco2, alpha, water_depth, facies)
    # TODO not used: I = precip .* facies.infiltration_coefficient #assume vertical infiltration
    reactive_surface =  facies.reactive_surface ./u"m^2/m^3"
    λ = precip .* facies.infiltration_coefficient ./ (alpha .* reactive_surface)
    eq = equilibrium(temp, pco2, precip, facies) # pass ceq Deq from the last function
    eq.denudation .* (1 - (λ ./ -water_depth) .* (1 - exp.(water_depth ./ λ))) * u"m/kyr"
end

function denudation(::Box{BT}, p::Dissolution, water_depth, slope, facies) where {BT<:Boundary}
    temp = p.temp ./ u"K"
    precip = p.precip ./u"m"
    pco2 = p.pco2 ./1.0u"atm"
    reactionrate = p.reactionrate ./u"m/yr"
    return (dissolution(temp, precip, pco2, reactionrate, water_depth, facies))
end

function redistribution(box::Box{BT}, p::Dissolution, water_depth, slope, facies) where {BT<:Boundary}
    return nothing
end

end
# ~/~ end