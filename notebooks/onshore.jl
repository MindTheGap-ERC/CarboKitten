### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ a2023b74-a36c-11ef-0520-b3ac38491c5a
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 02bf7fda-d6d4-4e13-ac30-e4f0c0f5821b
using Revise

# ╔═╡ b1348b07-5778-45e4-96be-be4f401b9431
using CarboKitten

# ╔═╡ 394b7989-ebd4-4aa2-8259-f52c701a2102
using CarboKitten.Models: ALCAP, OnshoreTransport as OT

# ╔═╡ 949e7a15-6c26-4cb5-85cb-249fdd24595a
using GLMakie

# ╔═╡ 26014b0c-b1a6-4c14-b996-5cc98af32a49
using CarboKitten.Visualization: summary_plot

# ╔═╡ fe8114d0-9776-4810-9f57-1a89e86c66dd
# ╠═╡ disabled = true
#=╠═╡
alcap_output = run_model(Model{ALCAP}, ALCAP.Example.INPUT, "../data/output/alcap.h5")
  ╠═╡ =#

# ╔═╡ d3353fee-3199-4a0c-baa7-bea0cb228b13
#=╠═╡
summary_plot(alcap_output)
  ╠═╡ =#

# ╔═╡ 5eb0ace5-c35f-4ffa-9081-3fd712960339


# ╔═╡ e1d2171e-eef4-423a-acf9-1cbb4de2b5bd
ot_output = let
	function ov(a, b)
		function (w)
			[(-max(0.0u"m/yr", b * (a - w) / a), 0.0u"m/yr"),
			 ((w <= a ? -b/a : 0.0u"yr^-1"), 0.0u"yr^-1")]
		end
	end
	
	FACIES = [
	    OT.Facies(
        maximum_growth_rate = 500u"m/Myr",
        extinction_coefficient = 0.8u"m^-1",
        saturation_intensity = 60u"W/m^2",
		diffusion_coefficient = 10u"m/yr",
		onshore_velocity = ov(10.0u"m", 0.0u"m/yr")),

	    OT.Facies(
        maximum_growth_rate = 400u"m/Myr",
        extinction_coefficient = 0.1u"m^-1",
        saturation_intensity = 60u"W/m^2",
		diffusion_coefficient = 50u"m/yr",
		onshore_velocity = ov(10.0u"m", -1.0u"m/yr")),

	    OT.Facies(
        maximum_growth_rate = 100u"m/Myr",
        extinction_coefficient = 0.005u"m^-1",
        saturation_intensity = 60u"W/m^2",
		diffusion_coefficient = 20u"m/yr",
		onshore_velocity = ov(10.0u"m", 1.0u"m/yr"))
	]

	function sea_level(t)
		10.0u"m" * sin(2π * t / 0.2u"Myr") + 3.0u"m" * sin(2π * t / 0.03u"Myr")
	end
	
	INPUT = OT.Input(
		tag = "ot1",
		box = CarboKitten.Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
		time = TimeProperties(
			Δt = 200.0u"yr",
			steps = 5000,
			write_interval = 10),
		sea_level = sea_level,
		initial_topography = (x, y) -> - x / 300.0,
		subsidence_rate = 50.0u"m/Myr",
		insolation = 400.0u"W/m^2",
		facies = FACIES,
		depositional_resolution = 0.5u"m",
		sediment_buffer_size = 50,
		disintegration_rate = 50.0u"m/Myr")

	run_model(Model{OT}, INPUT, "ot.h5")
end

# ╔═╡ 11fa7c68-40d0-49d8-b448-9db9626cc3dd
summary_plot(ot_output)

# ╔═╡ Cell order:
# ╠═a2023b74-a36c-11ef-0520-b3ac38491c5a
# ╠═02bf7fda-d6d4-4e13-ac30-e4f0c0f5821b
# ╠═b1348b07-5778-45e4-96be-be4f401b9431
# ╠═394b7989-ebd4-4aa2-8259-f52c701a2102
# ╠═949e7a15-6c26-4cb5-85cb-249fdd24595a
# ╠═26014b0c-b1a6-4c14-b996-5cc98af32a49
# ╠═fe8114d0-9776-4810-9f57-1a89e86c66dd
# ╠═d3353fee-3199-4a0c-baa7-bea0cb228b13
# ╠═5eb0ace5-c35f-4ffa-9081-3fd712960339
# ╠═e1d2171e-eef4-423a-acf9-1cbb4de2b5bd
# ╠═11fa7c68-40d0-49d8-b448-9db9626cc3dd
