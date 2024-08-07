

using Markdown
using InteractiveUtils

using HDF5

using Unitful

using GLMakie

using CarboKitten.Utility: in_units_of

using Statistics: mean

const Length = typeof(1.0u"m")

const Time = typeof(1.0u"Myr")

na = [CartesianIndex()]

struct Axes
	x::Vector{Length}
	y::Vector{Length}
	t::Vector{Time}
end

struct Header
	axes::Axes
	Δt::Time
	time_steps::Int
	height::Matrix{Float64}
	sea_level::Vector{Length}
	subsidence_rate::Float64
end

struct Data
	denudation::Array{Float64, 3}
	production::Array{Float64, 4}
	redistribution::Array{Float64, 3}
	sediment_elevation::Array{Float64, 3}
end

struct DataSlice
	denudation::Array{Float64, 2}
	production::Array{Float64, 3}
	redistribution::Array{Float64, 2}
	sediment_elevation::Array{Float64, 2}
end

function read_header(fid)
	attrs = HDF5.attributes(fid["input"])

	axes = Axes(
		fid["input/x"][]*u"m",
		fid["input/y"][]*u"m",
		fid["input/t"][]*u"Myr")
	
	return Header(
		axes,
		attrs["delta_t"][]*u"Myr",
		attrs["time_steps"][],
		fid["input/height"][]*u"m",
		fid["input/sea_level"][]*u"m",
		attrs["subsidence_rate"][]*u"m/Myr")
end

function read_data(filename)
	h5open(filename) do fid
		header = read_header(fid)
		data = Data(
			fid["denudation"][]*u"m",
			fid["production"][]*u"m",
			fid["redistribution"][]*u"m",
			fid["sediment_height"][]*u"m")
		header, data
	end
end

function read_slice(filename, slice...)
	h5open(filename) do fid
		header = read_header(fid)
		data = DataSlice(
			fid["denudation"][slice...]*u"m",
			fid["production"][slice...]*u"m",
			fid["redistribution"][slice...]*u"m",
			fid["sediment_height"][slice[2:end]...]*u"m")
		header, data
	end
end

read_slice("data/caps-test.h5", :, :, 25, :)

header, data = read_data("data/caps-test.h5")

η0 = header.height .- (header.axes.t[end] * header.subsidence_rate);

elevation(h::Header, d::Data) = let bl = h.height[:,:,na],
									sr = h.axes.t[end] * h.subsidence_rate
	cat(bl, bl .+ d.sediment_elevation; dims=3) .- sr
end

elevation(h::Header, d::DataSlice, y) = let bl = h.bedrock_elevation[:,y,na],
	       								 sr = h.axes.t[end] * h.subsidence_rate
	cat(bl, bl .+ d.sediment_elevation; dims=2) .- sr
end

η = elevation(header, data);

colormax(d::Data) = getindex.(argmax(d.deposition; dims=1)[1,:,:,:], 1)

colormax(d::DataSlice) = getindex.(argmax(d.deposition; dims=1)[1,:,:], 1)

colorrgb = eachslice(replace(data.deposition ./ sum(data.deposition; dims=1), NaN=>0.0), dims=(2,3,4)) .|> splat(RGBf);

function explode_quad_vertices(v::Array{Float64, 3})
	w, h, d = size(v)
	points = zeros(Float64, w, h-1, 2, d)
	n_vertices = 2 * w * (h-1)
	n_quads = (w - 1) * (h - 1)
	@views points[:, :, 1, :] = v[1:end, 1:end-1, :]
	@views points[:, :, 2, :] = v[1:end, 2:end, :]
	idx = reshape(1:n_vertices, w, (h-1), 2)
	vtx1 = reshape(idx[1:end-1, :, 1], n_quads)
	vtx2 = reshape(idx[2:end, :, 1], n_quads)
	vtx3 = reshape(idx[2:end, :, 2], n_quads)
	vtx4 = reshape(idx[1:end-1, :, 2], n_quads)
	return reshape(points, n_vertices, d), 
		vcat(hcat(vtx1, vtx2, vtx3), hcat(vtx1, vtx3, vtx4))
end

md"""
## Production profile
"""

let
	y = 25
	prod = sum(data.production; dims=1)[1,:,:,:]
	colormax = getindex.(argmax(data.production; dims=1)[1,:,:,:], 1)
	c = Array{Union{Int, Missing}}(missing, size(prod)...)
	mask = prod .!= 0.0u"m"
	c[mask] .= colormax[mask]
	heatmap(header.axes.x |> in_units_of(u"km"), header.axes.t[1:end-1] |> in_units_of(u"Myr"), c[:,y,:], colormap=:viridis)
end

md"""
## Detect production gaps
"""

local_height = η .- (header.subsidence_rate .* (header.axes.t .- header.axes.t[end]) .+ header.sea_level)[na, na, :];

let
	fig = Figure()
	ax = Axis(fig[1, 1])
	heatmap!(ax, local_height[:,25,:] .> 0u"m")
	fig
end

let
	x = header.axes.x |> in_units_of(u"km")
	t = header.axes.t |> in_units_of(u"Myr")
	lh = η .- (header.subsidence_rate .* (header.axes.t .- header.axes.t[end]) .+ header.sea_level)[na,na, :]
	fig = Figure(size=(800, 400))
	ax = Axis3(fig[1,1:3], azimuth=7π/4, zlabel="", xlabel="x(km)", ylabel="t(Myr)")
	surface!(ax, x, t, lh[:,25,:] |> in_units_of(u"m"))
	ax2 = Axis(fig[1,4:5], xlabel="t(Myr)", ylabel="-w(m)", title="x=1km")
	lines!(ax2, t, lh[8,25,:] |> in_units_of(u"m"))
	fig
end

function bean_counter(mask::BitArray{dim}) where {dim}
	visited = BitArray{dim}(undef, size(mask)...)
	visited .= false
	out = zeros(Int, size(mask)...)
	dxs = CartesianIndices(ntuple(_->3, dim)) .|> (x -> x - CartesianIndex(ntuple(_->2, dim))) |> filter(x->x!=CartesianIndex(ntuple(_->0, dim)...))
	group = 1
	
	for idx in CartesianIndices(mask)
		visited[idx] && continue
		visited[idx] = true
		mask[idx] || continue
		out[idx] = group

		stack = idx .+ dxs
		while !isempty(stack)
			jdx = pop!(stack)
			checkbounds(Bool, mask, jdx) || continue
			visited[jdx] && continue
			visited[jdx] = true
			mask[jdx] || continue
			out[jdx] = group
			append!(stack, jdx .+ dxs)
		end
		group += 1
	end
	return out, group-1
end

bean_counter(BitArray([0 1 0; 0 0 0; 1 1 0]))

gaps, n_gaps = bean_counter(local_height .> 0u"m")

heatmap(gaps[:, 25, :] .|> (x-> x==0 ? missing : x))

md"""
## Sediment profile
"""

function plot_sediment_profile(filename, y)
	header, data = read_slice(filename, :, :, y, :)
	x = header.axes.x |> in_units_of(u"km")
	t = header.axes.t |> in_units_of(u"Myr")
	ξ = elevation(header, data, y)  # |> in_units_of(u"m")
	
	verts = zeros(Float64, length(x), length(t), 2)
	@views verts[:, :, 1] .= x
	@views verts[:, :, 2] .= ξ |> in_units_of(u"m")
	v, f = explode_quad_vertices(verts)

	water_depth = ξ .- (header.subsidence_rate .* (header.axes.t .- header.axes.t[end]) .+ header.sea_level)[na, :]
	gaps, n_gaps = bean_counter(water_depth .> 0u"m")

	fig = Figure(size=(800,600))
	ax = Axis(fig[1,1])
	c = reshape(colormax(data)[:, :], length(x) * (length(t) - 1))
	mesh!(ax, v, f, color=vcat(c, c), alpha=1.0)

	for g = 1:n_gaps
		size = sum(gaps .== g)
		if size < 1000
			continue
		end
		gap = mean.(skipmissing.(eachslice(CartesianIndices(ξ) .|> (i -> gaps[i] == g ? ξ[i] : missing), dims=(1,))))
		lines!(ax, x, gap |> in_units_of(u"m"), color=:white, linewidth=2, linestyle=:dash)
	end
	
	fig
end

plot_sediment_profile("../data/test.h5", 25)

md"""
## Sediment profile with time

Seen from above this is a Wheeler diagram, from the side we get a sediment profile.
"""

let
	y = 25
	x = header.axes.x |> in_units_of(u"km")
	t = header.axes.t |> in_units_of(u"Myr")

	verts = zeros(Float64, length(x), length(t), 3)
	@views verts[:, :, 1] .= (header.axes.x |> in_units_of(u"km"))
	@views verts[:, :, 3] .= (η[:, y, :] |> in_units_of(u"m"))
	@views verts[:, :, 2] .= (header.axes.t |> in_units_of(u"Myr"))[na,:]
	v, f = explode_quad_vertices(verts)
	c = reshape(colormax(data)[:, y, :], length(x) * (length(t) - 1))

	fig = Figure()
	ax = Axis3(fig[1,1], azimuth=7π/4, xlabel="x(km)", ylabel="t(Myr)", zlabel="z(m)")
	mesh!(ax, v, f, color=vcat(c, c), alpha=1.0)
	save("active-layer-test.png", fig)
	fig
end

md"""
## Top layer
"""

let
	fig = Figure()
	ax = Axis3(fig[1, 1], azimuth=7π/4)
	surface!(ax, header.axes.x |> in_units_of(u"km"), header.axes.y |> in_units_of(u"km"), η[:,:,end] |> in_units_of(u"m"), color=colorrgb[:,:,end])
	fig
end

