# ~/~ begin <<docs/src/ca-with-production.md#examples/caps-osc.jl>>[init]
using CarboKitten.CaProd

DEFAULT_INPUT = CaProd.Input(
    sea_level = t -> 4 * sin(2π * t / 0.2), 
    subsidence_rate = 50.0,
    initial_depth = x -> x / 2,
    grid_size = (50, 100),
    phys_scale = 1.0,
    Δt = 0.0001,
    write_interval = 10,
    time_steps = 10000,
    facies = [
        CaProd.Facies((4, 10), (6, 10), 500.0, 0.8, 300, 1000, 2730, 0.5),
        CaProd.Facies((4, 10), (6, 10), 400.0, 0.1, 300, 1000, 2730, 0.5),
        CaProd.Facies((4, 10), (6, 10), 100.0, 0.005, 300, 1000, 2730, 0.5)
    ],
    insolation = 2000.0,
    temp = 288.0,
    precip = 1000.0,
    pco2 = 10^(-1.5),
    alpha = 2e-6

)

CaProd.main(DEFAULT_INPUT, "data/caps-osc.h5")
# ~/~ end