# Denudation component

``` {.julia file=src/Components/Denudation.jl}
@compose module Denudation
    @mixin Boxes, WaterDepth, SedimentBuffer, FaciesBase
    using ..Common
    using CarboKitten.Denudation: DenudationType
    using CarboKitten.BoundaryTrait
    using CarboKitten.Stencil: stencil
    using CarboKitten.Denudation: Dissolution, EmpiricalDenudation, PhysicalErosion, NoDenudation, denudation, redistribution, slope_kernel

    export slope_function
    export Dissolution, EmpiricalDenudation, PhysicalErosion, NoDenudation, denudation, redistribution, slope_kernel


    @kwdef struct Facies <: AbstractFacies
        reactive_surface::typeof(1.0u"m^2/m^3") #reactive surface
        mass_density::typeof(1.0u"kg/m^3") #density of different carb factory
        infiltration_coefficient::Float64 #infiltration coeff
        erodibility::typeof(1.0u"m/yr")
    end

    @kwdef struct Input <: AbstractInput
        denudation::DenudationType
    end

    function slope_function(input,box::Box{BT}) where {BT<:Boundary}
        return stencil(Float64, BT, (3, 3), slope_kernel) 
    end
end
```
