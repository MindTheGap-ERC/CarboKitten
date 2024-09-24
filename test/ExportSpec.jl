# ~/~ begin <<docs/src/data-export.md#test/ExportSpec.jl>>[init]
using CarboKitten.Export: Axes, Header, Data, data_export, CSVExportTrait
using CSV: read as read_csv
using DataFrames
using Unitful

const Amount = typeof(1.0u"m")

const AXES1 = Axes(
  x = [0.0, 1.0, 2.0] * u"m",
  y = [1.0] * u"m",
  t = (0.0:0.1:1.0) * u"Myr")

const HEADER1 = Header(
  axes = AXES1,
  Δt = 0.1u"Myr",
  time_steps = 10,
  bedrock_elevation = zeros(typeof(1.0u"m"), 3, 3),
  sea_level = zeros(typeof(1.0u"m"), 10),
  subsidence_rate = 10u"m/Myr")

# 1x3x1x10
# case1: stable production, no disintegration
# case2: stable production, top-hat disintegration
# case3: linear production, no disintegration

const PRODUCTION1 = reshape(
  hcat(ones(Amount, 10),
       ones(Amount, 10),
       cumsum(ones(Amount, 10)) / 5.5)',
  1, 3, 1, 10)

const DISINTEGRATION1 = reshape(
  hcat(zeros(Amount, 10),
       1:10 .|> (x -> x < 4 || x > 6 ? 0.0u"m" : 2.0u"m"),
       zeros(Amount, 10))',
  1, 3, 1, 10)

const ELEVATION1 = cumsum(PRODUCTION1 .- DISINTEGRATION1; dims=4)[1,:,:,:]

const DATA1 = Data(
  disintegration = DISINTEGRATION1,
  production = PRODUCTION1,
  deposition = PRODUCTION1 .- DISINTEGRATION1,
  sediment_elevation = ELEVATION1)

const GRID_LOCATIONS1 = [(1,1), (2, 1), (3,1)]

@testset "Data Export" begin
  @testset "Hither and Dither" begin
      buffer = UInt8[]
      io = IOBuffer(buffer, read=true, write=true)
      data_export(CSVExportTrait{:sediment_accumulation_curve}, io, HEADER1, DATA1, GRID_LOCATIONS1)
      println(read(IOBuffer(buffer), String))
      df = read_csv(IOBuffer(buffer), DataFrame)
      rename!(df, (n => split(n)[1] for n in names(df))...)
      @test df.sac1 ≈ ELEVATION1[1, 1, :] / u"m"
      @test df.sac2 ≈ ELEVATION1[2, 1, :] / u"m"
      @test df.sac3 ≈ ELEVATION1[3, 1, :] / u"m"
  end
end
# ~/~ end
