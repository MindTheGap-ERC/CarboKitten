# ~/~ begin <<docs/src/model-alcap.md#examples/model/alcap/run.jl>>[init]
#| requires: src/Model/ALCAP2.jl
#| creates: data/output/alcap2.h5

module Script

using Unitful
using CarboKitten.Components
using CarboKitten.Components.Common
using CarboKitten.Components.DenudationConfig
using CarboKitten.Model: ALCAP2 as ALCAP
using CarboKitten.Model: WithDenudation2 as WDn
using CarboKitten.Export: data_export, CSV
using CarboKitten.Denudation

const m = u"m"
const Myr = u"Myr"

const PATH = "data/output"
const TAG = "empirical"

const FACIES = [
    WDn.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=10000u"m",
        reactive_surface=1000u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.23u"m/yr"
        ),
    WDn.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=5000u"m",
        reactive_surface=1000u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.23u"m/yr"
        ),
    WDn.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=7000u"m",
        reactive_surface=1000u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.23u"m/yr"
        )
]

const PERIOD = 0.2Myr
const AMPLITUDE = 20.0m
const DENUDATION = EmpiricalDenudation(precip=1.0u"m")


const INPUT = WDn.Input(
    tag="$TAG",
    box=Box{Shelf}(grid_size=(100, 50), phys_scale=150.0m),
    time=TimeProperties(
        Δt=0.0002Myr,
        steps=5000,
        write_interval=1),
    ca_interval=1,
    bedrock_elevation=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate=20.0m / Myr,
    disintegration_rate=500.0m / Myr,
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5m,
    facies=FACIES,
    denudation = DENUDATION)

function main()
    H5Writer.run(Model{WDn}, INPUT, "$(PATH)/$(TAG).h5")

    data_export(
        CSV(tuple.(10:20:70, 25),
          :sediment_accumulation_curve => "$(PATH)/$(TAG)_sac.csv",
          :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
          :stratigraphic_column => "$(PATH)/$(TAG)_sc.csv",
          :metadata => "$(PATH)/$(TAG).toml"),
        "$(PATH)/$(TAG).h5")
end

end

Script.main()
# ~/~ end