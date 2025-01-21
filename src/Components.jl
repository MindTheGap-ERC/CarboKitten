# ~/~ begin <<docs/src/components/components.md#src/Components.jl>>[init]
module Components

export Tag, TimeIntegration, Boxes, WaterDepth, FaciesBase, Production,
    CAProduction, CellularAutomaton, H5Writer, ActiveLayer, SedimentBuffer, Denudation

using ModuleMixins: @compose

include("Components/Common.jl")
include("Components/Tag.jl")
include("Components/TimeIntegration.jl")
include("Components/Boxes.jl")
include("Components/WaterDepth.jl")
include("Components/FaciesBase.jl")
include("Components/Production.jl")
include("Components/CellularAutomaton.jl")
include("Components/CAProduction.jl")

include("Components/SedimentBuffer.jl")
include("Components/ActiveLayer.jl")
include("Components/Denudation.jl")

include("Components/H5Writer.jl")

list_components() = filter(
    c->:AST in names(c, all=true),
    map(eval, names(Components)))

end
# ~/~ end
