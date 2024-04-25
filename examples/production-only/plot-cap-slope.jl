# ~/~ begin <<docs/src/ca-with-production.md#examples/production-only/plot-cap-slope.jl>>[init]
#| creates: docs/src/fig/b13-crosssection.png
#| requires: data/ca-prod-slope.h5 ext/VisualizationExt.jl
#| collect: figures

module Script
# using CarboKitten.Visualization
using CairoMakie
using Unitful
using HDF5
using GeometryBasics

module Data

using Unitful
using HDF5
using AxisArrays

@kwdef struct CKData
    Δt::typeof(1.0u"Myr")
    subsidence_rate::typeof(1.0u"m/Myr")
    t_axis::typeof(Float64[]u"Myr")
    x_axis::typeof(Float64[]u"m")
    initial_height::typeof(Float64[]u"m")
    sediment::typeof(AxisArray(Array{Float64,4}(undef,0,0,0,0)u"m/Myr", Axis{:x}(), Axis{:y}(), Axis{:facies}(), Axis{:t}()))
end


function read_data(datafile)
    h5open(datafile, "r") do fid
        attr = HDF5.attributes(fid["input"])
        t_axis = fid["input/t"][]u"Myr"
        x_axis = fid["input/x"][]u"m"
        y_axis = fid["input/y"][]u"m"
        CKData(
            Δt = attr["delta_t"][]u"Myr",
            subsidence_rate = attr["subsidence_rate"][]u"m/Myr",
            t_axis = fid["input/t"][]u"Myr",
            x_axis = fid["input/x"][]u"m",
            initial_height = fid["input/height"][]u"m",
            sediment = AxisArray(fid["sediment"][]u"m/Myr", Axis{:x}(x_axis), Axis{:y}(y_axis), Axis{:facies}(1:n_facies), Axis{:t}(t_axis))
        )
    end
end

end # module Data

function plot_crosssection(pos, datafile)
    # x: 1-d array with x-coordinates
    # t: 1-d array with time-coordinates (n_steps + 1)
    # h[x, t]: height fn, monotonic increasing in time
    # p[x, facies, t]: production rate
    # taken at y = y_max / 2, h[x, 1] is initial height
    n_facies, x, t, h, p = h5open(datafile, "r") do fid
        attr = HDF5.attributes(fid["input"])
        Δt = attr["delta_t"][]
        subsidence_rate = attr["subsidence_rate"][]
        t_end = fid["input/t"][end-1]
        total_subsidence = subsidence_rate * t_end
        total_sediment = sum(fid["sediment"][]; dims=3)
        initial_height = fid["input/height"][]
        center = div(size(total_sediment)[1], 2)
        elevation = cumsum(total_sediment; dims=4)[:, center, 1, :] .* Δt .- initial_height .- total_subsidence
        t = fid["input/t"][]
        n_facies = size(fid["sediment"])[3]

        return n_facies,
        fid["input/x"][],
        [t; Δt * attr["time_steps"][]],
        hcat(.-initial_height .- total_subsidence, elevation),
        fid["sediment"][:, center, :, :]
    end

    pts = vec(Point{2,Float64}.(x, h[:, 2:end]))
    c = vec(argmax(p; dims=2)[:, 1, :] .|> (c -> c[2]))
    rect = Rect2(0.0, 0.0, 1.0, 1.0)
    m_tmp = GeometryBasics.mesh(Tesselation(rect, (100, 1000)))
    m = GeometryBasics.Mesh(pts, faces(m_tmp))

    # pts = vec(Point{2,Float64}.(x, h))
    # c = argmax(p; dims=2)[:,1,:] .|> (c -> c[2])
    # w = size(x)[1]

    # face(idx) = let k = idx[1] + idx[2]*w
    #     TriangleFace(k, k+1, k+1+w), TriangleFace(k+1+w, k+w, k)
    # end

    ax = Axis(pos, xlabel="location", ylabel="depth", limits=((-12, x[end]), nothing))
    # for f in 1:n_facies
    #     locs = CartesianIndices((size(x)[1], size(t)[1] - 1))[c .== f]
    #     triangles = collect(Iterators.flatten(face.(locs)))
    #     m = GeometryBasics.Mesh(pts, triangles)
    #     mesh!(ax, m)
    # end

    mesh!(ax, m, color=c, alpha=0.7)
    for idx in [1, 501, 1001]
        lines!(ax, x, h[:, idx], color=:black)
        text!(ax, -2.0, h[1, idx]; text="$(t[idx]) Myr", align=(:right, :center))
    end
    for idx in [250, 750]
        lines!(ax, x, h[:, idx], color=:black, linewidth=0.5)
    end
end

function main()
    f = Figure()
    plot_crosssection(f[1, 1], "data/ca-prod-slope.h5")
    save("docs/src/fig/b13-crosssection.png", f)
end
end

# Script.main()
# ~/~ end