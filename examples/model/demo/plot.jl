# ~/~ begin <<docs/src/model-alcap.md#examples/model/alcap/plot.jl>>[init]
#| creates: docs/src/_fig/alcaps-alternative.png
#| requires: data/output/alcap2.h5
#| collect: figures

using GLMakie
using CarboKitten.Visualization

save("docs/src/_fig/alcaps-la04.png", summary_plot("data/output/alcap2.h5"))
# ~/~ end