# ~/~ begin <<docs/src/data-export.md#test/ExportSpec.jl>>[init]
using CarboKitten.Export: Axes, Header, Data, data_export, CSVExportTrait,
    age_depth_model, extract_sac, extract_sc, CSV
using CSV: read as read_csv
using TOML
using DataFrames
using Unitful

const Amount = typeof(1.0u"m")

# ~/~ begin <<docs/src/data-export.md#export-test-case>>[init]
const AXES1 = Axes(
    x=[0.0, 1.0, 2.0] * u"m",
    y=[1.0] * u"m",
    t=(0.0:0.1:1.0) * u"Myr")

const HEADER1 = Header(
    tag="test",
    axes=AXES1,
    write_interval=1,
    Δt=0.1u"Myr",
    time_steps=10,
    initial_topography=zeros(typeof(1.0u"m"), 3, 3),
    sea_level=zeros(typeof(1.0u"m"), 10),
    subsidence_rate=10u"m/Myr")

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

const ELEVATION1 = cat(
    [0.0, 0.0, 0.0]u"m",
    cumsum(PRODUCTION1 .- DISINTEGRATION1; dims=4)[1, :, :, :];
    dims=3)

const DATA1 = Data(
    disintegration=DISINTEGRATION1,
    production=PRODUCTION1,
    deposition=PRODUCTION1 .- DISINTEGRATION1,
    sediment_elevation=ELEVATION1)

const GRID_LOCATIONS1 = [(1, 1), (2, 1), (3, 1)]
# ~/~ end

@testset "Data Export" begin
    # ~/~ begin <<docs/src/data-export.md#export-test>>[init]
    @testset "Hither and Dither" begin
        io = IOBuffer(UInt8[], read=true, write=true)
        data_export(CSVExportTrait{:sediment_accumulation_curve}, io, HEADER1, DATA1, GRID_LOCATIONS1)
        seek(io, 0)
        df = read_csv(io, DataFrame)
        rename!(df, (n => split(n)[1] for n in names(df))...)
        @test df.sac1 ≈ ELEVATION1[1, 1, :] / u"m"
        @test df.sac2 ≈ ELEVATION1[2, 1, :] / u"m"
        @test df.sac3 ≈ ELEVATION1[3, 1, :] / u"m"
    end
    # ~/~ end
    # ~/~ begin <<docs/src/data-export.md#export-test>>[1]
    @testset "ADM Monotonicity" begin
        sac = extract_sac(HEADER1, DATA1, GRID_LOCATIONS1)
        adm = sac |> age_depth_model

        @test sac.sac1 == adm.adm1
        @test sac.sac3 == adm.adm3
        @test sac.sac2 != adm.adm2

        @test all(adm.adm2[2:end] .- adm.adm2[1:end-1] .>= 0.0u"m")
    end
    # ~/~ end
    # ~/~ begin <<docs/src/data-export.md#export-test>>[2]
    @testset "SC sum equals ADM" begin
        sac = extract_sac(HEADER1, DATA1, GRID_LOCATIONS1)
        adm = sac |> age_depth_model
        sc = extract_sc(HEADER1, DATA1, GRID_LOCATIONS1)
        @test [0.0u"m"; cumsum(sc.sc1_f1)] ≈ adm.adm1
        @test [0.0u"m"; cumsum(sc.sc2_f1)] ≈ adm.adm2
        @test [0.0u"m"; cumsum(sc.sc3_f1)] ≈ adm.adm3
    end
    # ~/~ end

    @testset "Write to folder" begin
        mktempdir() do path
            spec = CSV(GRID_LOCATIONS1,
                :sediment_accumulation_curve => joinpath(path, "sac.csv"),
                :age_depth_model => joinpath(path, "adm.csv"),
                :stratigraphic_column => joinpath(path, "sc.csv"),
                :water_depth => joinpath(path, "wd.csv"),
                :metadata => joinpath(path, "metadata.toml"))
            data_export(spec, HEADER1, DATA1)
            for f in values(spec.output_files)
                @test isfile(f)
            end

            metadata = TOML.parsefile(spec.output_files[:metadata])
            @test IdDict(Symbol(k) => v for (k, v) in metadata["files"]) == spec.output_files
            @test length(metadata["locations"]) == 3
            adm = read_csv(spec.output_files[:age_depth_model], DataFrame)
            rename!(adm, (n => split(n)[1] for n in names(adm))...)
            @test adm == ustrip(extract_sac(HEADER1, DATA1, GRID_LOCATIONS1) |> age_depth_model)
        end
    end
end
# ~/~ end
