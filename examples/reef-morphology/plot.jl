# ~/~ begin <<docs/src/model-alcap.md#examples/model/alcap/plot.jl>>[init]
#| creates: docs/src/_fig/alcaps-alternative.png
#| requires: data/output/alcap-example.h5
#| collect: figures

using GLMakie
using CarboKitten.Visualization

save("docs/src/_fig/alcaps-ramp.png", summary_plot("data/output/alcap-ramp.h5"))
# ~/~ end
