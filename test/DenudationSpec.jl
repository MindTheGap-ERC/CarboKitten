using Test
using Unitful

@testset "DenudationTST" begin
    import CarboKitten.Denudation.DissolutionMod: dissolution
    import CarboKitten.Denudation.EmpiricalDenudationMod: empirical_denudation, slope_kernel
    import CarboKitten.Denudation.PhysicalErosionMod: physical_erosion, mass_erosion, total_mass_redistribution
    using CarboKitten.Stencil: Periodic, Reflected, stencil
    using CarboKitten.Config: Vectors, TimeProperties
    using CarboKitten.Boxes: Box
    using CarboKitten.Burgess2013.CA: step_ca, run_ca
    using CarboKitten.Models.WithDenudation: Input, Facies
    using CarboKitten.Denudation: denudation, redistribution, Dissolution, NoDenudation, PhysicalErosion, EmpiricalDenudation


    DENUDATION_HIGH_CO2 = Dissolution(temp = 293.0u"K",precip = 1.0u"m", pco2 = 10^(-1.5)*u"atm",reactionrate = 2e-3u"m/yr")
    DENUDATION_LOW_CO2 = Dissolution(temp = 293.0u"K",precip = 1.0u"m", pco2 = 10^(-2.5)*u"atm",reactionrate = 2e-3u"m/yr")
    DENUDATION_LOW_P = EmpiricalDenudation(precip = 0.8u"m")
    DENUDATION_HIGH_P = EmpiricalDenudation(precip = 1.0u"m")
    DENUDATION_PHYS = PhysicalErosion()
    MODEL1 = [
        Facies(viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate = 500u"m/Myr",
        extinction_coefficient = 0.8u"m^-1",
        saturation_intensity = 60u"W/m^2",
        diffusion_coefficient=10000u"m",
        reactive_surface = 1000u"m^2/m^3",
        mass_density = 2730u"kg/m^3",
        infiltration_coefficient= 0.5,
        erodibility = 0.23u"m/yr"),

        Facies(viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate = 400u"m/Myr",
        extinction_coefficient = 0.1u"m^-1",
        saturation_intensity = 60u"W/m^2",
        diffusion_coefficient=10000u"m",
        reactive_surface = 1000u"m^2/m^3",
        mass_density = 2730u"kg/m^3",
        infiltration_coefficient= 0.5,
        erodibility = 0.23u"m/yr"),

        Facies(viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate = 100u"m/Myr",
        extinction_coefficient = 0.005u"m^-1",
        saturation_intensity = 60u"W/m^2",
        diffusion_coefficient=10000u"m",
        reactive_surface = 1000u"m^2/m^3",
        mass_density = 2730u"kg/m^3",
        infiltration_coefficient= 0.5,
        erodibility = 0.23u"m/yr")
    ]

    box = Box{Periodic{2}}(grid_size=(5, 5), phys_scale=1.0u"km")
    n_facies = length(MODEL1)
    ca_init = [ 0  0  1  3  3
                0  1  3  2  1
                2  0  1  0  1
                1  3  3  3  0
                1  3  2  3  2]

    struct test_state
        ca::Array{Int}
    end

    STATE1 = test_state(ca_init)

    denudation_mass_HIGH_CO2 = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    denudation_mass_LOW_CO2 = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    denudation_mass_LOW_P = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    denudation_mass_HIGH_P = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    denudation_mass_phys = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    denudation_mass_phys_flat = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    redistribution_mass = zeros(typeof(0.0u"m"),n_facies,box.grid_size...)

    water_depth = -100 .* [ 0.989943  0.48076   0.518983  0.997996   0.895681
                    0.872733  0.208779  0.882917  0.550494   0.674066
                    0.57987   0.619433  0.769506  0.593786   0.856186
                    0.407728  0.469545  0.896348  0.473817   0.797112
                    0.610194  0.921632  0.322729  0.0103646  0.691191]
    water_depth_flat = -0.5 .* ones(box.grid_size...)
    inf_map = ones(box.grid_size...)
    slope = rand(Float64, box.grid_size...)
    slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel)
    slopefn(water_depth, slope, box.phys_scale ./u"m")
    slope_flat = zeros(box.grid_size...)
    slopefn(water_depth_flat, slope_flat, box.phys_scale ./u"m")

    (denudation_mass_HIGH_CO2) = denudation(box, DENUDATION_HIGH_CO2, water_depth, slope,MODEL1,STATE1)
    (denudation_mass_LOW_CO2) = denudation(box, DENUDATION_LOW_CO2, water_depth, slope,MODEL1,STATE1)
    (denudation_mass_LOW_P) = denudation(box, DENUDATION_LOW_P, water_depth, slope,MODEL1,STATE1)
    (denudation_mass_HIGH_P) = denudation(box, DENUDATION_HIGH_P, water_depth, slope,MODEL1,STATE1)

    (denudation_mass_phys) = denudation(box, DENUDATION_PHYS, water_depth, slope,MODEL1,STATE1)
    (denudation_mass_phys_flat) = denudation(box, DENUDATION_PHYS, water_depth_flat, slope_flat,MODEL1,STATE1)
    for idx in CartesianIndices(STATE1.ca)
        f = STATE1.ca[idx]
        if f == 0
            continue
        end

    inf_map[idx] = MODEL1[f].infiltration_coefficient
    end

    (redistribution_mass) = redistribution(box,DENUDATION_PHYS,denudation_mass_phys .*1.0u"Myr",water_depth)
    # println(redistribution_mass)
    # println(denudation_mass_phys.*1.0u"Myr")
    @test sum(denudation_mass_HIGH_CO2) > sum(denudation_mass_LOW_CO2)
    @test sum(denudation_mass_LOW_P) < sum(denudation_mass_HIGH_P)
    @test sum(denudation_mass_phys) > sum(denudation_mass_phys_flat)
    @test sum(denudation_mass_phys .*1.0u"Myr") ≈ sum(redistribution_mass)

    #regression_test
    #@test 1.8*0.95 < abs.(sum(denudation_mass_HIGH_CO2)) ./u"m/kyr" < 1.8 *1.05
    #@test 477*0.95 < abs.(sum(denudation_mass_LOW_P)) ./u"m/kyr" < 477*1.05
    #@test 0.15* 0.7115 < sum(denudation_mass_phys) ./u"m/kyr" < 0.15*1.05

end
