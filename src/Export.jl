# ~/~ begin <<docs/src/data-export.md#src/Export.jl>>[init]
module Export

export Data, DataSlice, DataColumn, Header, CSV, read_data, read_slice, read_column, data_export

using HDF5
import CSV: write as write_csv
using TOML

using Unitful
using DataFrames
using .Iterators: flatten

const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")
const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")

const na = [CartesianIndex()]

# ~/~ begin <<docs/src/data-export.md#export-specification>>[init]
abstract type ExportSpecification end

@kwdef struct CSV <: ExportSpecification
    grid_locations::Vector{NTuple{2,Int}}
    output_files::IdDict{Symbol,String}
end

CSV(grid_locations, kwargs...) = CSV(grid_locations, IdDict(kwargs...))
# ~/~ end

@kwdef struct Axes
    x::Vector{Length}
    y::Vector{Length}
    t::Vector{Time}
end

@kwdef struct Header
    tag::String
    axes::Axes
    Δt::Time
    write_interval::Int
    time_steps::Int
    initial_topography::Matrix{Amount}
    sea_level::Vector{Length}
    subsidence_rate::Rate
end

@kwdef struct Data
    disintegration::Array{Amount,4}
    production::Array{Amount,4}
    deposition::Array{Amount,4}
    sediment_elevation::Array{Amount,3}
end

struct DataSlice
    slice::NTuple{2,Union{Colon,Int}}
    disintegration::Array{Amount,3}
    production::Array{Amount,3}
    deposition::Array{Amount,3}
    sediment_elevation::Array{Amount,2}
end

struct DataColumn
    slice::NTuple{2,Int}
    disintegration::Array{Amount,2}
    production::Array{Amount,2}
    deposition::Array{Amount,2}
    sediment_elevation::Array{Amount,1}
end

function read_header(fid)
    attrs = HDF5.attributes(fid["input"])

    axes = Axes(
        fid["input/x"][] * u"m",
        fid["input/y"][] * u"m",
        fid["input/t"][] * u"Myr")

    return Header(
        attrs["tag"][],
        axes,
        attrs["delta_t"][] * u"Myr",
        attrs["write_interval"][],
        attrs["time_steps"][],
        fid["input/initial_topography"][] * u"m",
        fid["input/sea_level"][] * u"m",
        attrs["subsidence_rate"][] * u"m/Myr")
end

function read_data(filename)
    h5open(filename) do fid
        header = read_header(fid)
        data = Data(
            fid["disintegration"][] * u"m",
            fid["production"][] * u"m",
            fid["deposition"][] * u"m",
            fid["sediment_height"][] * u"m")
        header, data
    end
end

read_slice(fid::HDF5.File, slice...) = DataSlice(
    slice,
    fid["disintegration"][:, slice..., :] * u"m",
    fid["production"][:, slice..., :] * u"m",
    fid["deposition"][:, slice..., :] * u"m",
    fid["sediment_height"][slice..., :] * u"m")

function read_slice(filename::AbstractString, slice...)
    h5open(filename) do fid
        header = read_header(fid)
        data = read_slice(fid, slice...)
        header, data
    end
end

read_column(fid::HDF5.File, slice...) = DataColumn(
    slice,
    fid["disintegration"][:, slice..., :] * u"m",
    fid["production"][:, slice..., :] * u"m",
    fid["deposition"][:, slice..., :] * u"m",
    fid["sediment_height"][slice..., :] * u"m")

function read_column(filename::AbstractString, slice...)
    h5open(filename) do fid
        header = read_header(fid)
        data = read_column(fid, slice...)
        header, data
    end
end

# ~/~ begin <<docs/src/data-export.md#export-function>>[init]
"""
    unitful_headers(df::DataFrame)

Gets a string representation for all column names including their unit.
Returns a `Vector{String}`.
"""
unitful_headers(df::DataFrame) =
    ["$(e.variable) [$(unit(e.eltype))]" for e in eachrow(describe(df))]

"""
    ustrip(df::DataFrame)

Strip units from a `DataFrame`. Returns a new `DataFrame`.
"""
Unitful.ustrip(df::DataFrame) =
    DataFrame((e.variable => df[!, e.variable] / unit(e.eltype)
               for e in eachrow(describe(df)))...)

"""
    write_unitful_csv(io::IO, df::DataFrame)

Write a CSV from a `DataFrame` with `Unitful` units. The units will be
represented in the CSV header, and stripped from the individual values.
"""
write_unitful_csv(io, df::DataFrame) =
    write_csv(io, ustrip(df), header=unitful_headers(df))
# ~/~ end
# ~/~ begin <<docs/src/data-export.md#export-function>>[1]
Base.accumulate(f) = (args...; kwargs...) -> accumulate(f, args...; kwargs...)

"""
    age_depth_model(sediment_accumulation_curve::Vector)
    age_depth_model(sediment_accumulation_curve::DataFrame)

Compute the ADM from the SAC. Implemented as:

    reverse ∘ accumulate(min) ∘ reverse

The `DataFrame` version `select`s SAC columns, transformed into ADM.
"""
age_depth_model(sac::Vector{T}) where {T} = sac |> reverse |> accumulate(min) |> reverse
age_depth_model(sac_df::DataFrame) =
    let sac_cols = filter(contains("sac"), names(sac_df)),
        adm_cols = replace.(sac_cols, "sac" => "adm")

        select(sac_df, "time", (sac => age_depth_model => adm
                                for (sac, adm) in zip(sac_cols, adm_cols))...)
    end
# ~/~ end
# ~/~ begin <<docs/src/data-export.md#export-function>>[2]
"""
    stratigraphic_column(header::Header, data::Data, loc::NTuple{2,Int}, facies::Int)

Compute the Stratigraphic Column for a given grid position `loc` and `facies` index.
Returns an `Array{Quantity, 2}` where the `Quantity` is in units of meters.
"""
function stratigraphic_column(header::Header, data::Data, loc::NTuple{2,Int}, facies::Int)
    dc = DataColumn(
        loc,
        data.disintegration[:, loc..., :],
        data.production[:, loc..., :],
        data.deposition[:, loc..., :],
        data.sediment_elevation[loc..., :])
    return stratigraphic_column(header, dc, facies)
end

function stratigraphic_column(header::Header, data::DataColumn, facies::Int)
    n_times = length(header.axes.t) - 1
    sc = zeros(typeof(1.0u"m"), n_times)

    for ts = 1:n_times
        acc = data.deposition[facies, ts] - data.disintegration[facies, ts]
        if acc > 0.0u"m"
            sc[ts] = acc
            continue
        end
        ts_down = ts - 1
        while acc < 0.0u"m"
            ts_down < 1 && break
            if -acc < sc[ts_down]
                sc[ts_down] -= acc
                break
            else
                acc += sc[ts_down]
                sc[ts_down] = 0.0u"m"
            end
            ts_down -= 1
        end
    end

    sc
end
# ~/~ end
# ~/~ begin <<docs/src/data-export.md#export-function>>[3]
struct CSVExportTrait{S} end

function data_export(spec::T, filepath::String) where {T<:ExportSpecification}
    data_export(spec, read_data(filepath)...)
end

function data_export(spec::CSV, header::Header, data::Data)
    for (key, filename) in spec.output_files
        if key == :metadata
            md = Dict(
                "global" => Dict(
                    "tag" => header.tag,
                    "subsidence_rate" => header.subsidence_rate,
                    "time_steps" => header.time_steps,
                    "delta_t" => header.Δt),
                "locations" => [Dict(
                    "number" => i,
                    "x" => header.axes.x[loc[1]],
                    "y" => header.axes.y[loc[2]],
                    "initial_topography" => header.initial_topography[loc...])
                                for (i, loc) in enumerate(spec.grid_locations)],
                "files" => spec.output_files)
            open(filename, "w") do io
                TOML.print(io, md) do obj
                    if obj isa Quantity
                        [ustrip(obj), string(unit(obj))]
                    else
                        obj
                    end
                end
            end
            continue
        end
        open(filename, "w") do io
            data_export(CSVExportTrait{key}, io, header, data, spec.grid_locations)
        end
    end
end

function data_export(::Type{CSVExportTrait{S}}, args...) where {S}
    error("Unknown CSV data export: `$(S)`")
end

function data_export(::Type{CSVExportTrait{:sediment_accumulation_curve}},
    io::IO, header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})

    sac = extract_sac(header, data, grid_locations)
    write_unitful_csv(io, sac)
end

function data_export(::Type{CSVExportTrait{:age_depth_model}},
    io::IO, header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})

    adm = extract_sac(header, data, grid_locations) |> age_depth_model
    write_unitful_csv(io, adm)
end

function data_export(::Type{CSVExportTrait{:stratigraphic_column}},
    io::IO, header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})

    sc = extract_sc(header, data, grid_locations)
    write_unitful_csv(io, sc)
end

function data_export(::Type{CSVExportTrait{:water_depth}},
    io::IO, header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})
    wd = extract_wd(header, data, grid_locations)
    write_unitful_csv(io, wd)
end

"""
    extract_sac(header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})

Extract Sediment Accumumlation Curve (SAC) from the data. The SAC is directly copied from
`data.sediment_elevation`. Returns a `DataFrame` with `time` and `sac<n>` columns where `<n>`
is in the range `1:length(grid_locations)`.
"""
function extract_sac(header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})
    DataFrame(:time => header.axes.t[1:end],
        (Symbol("sac$(i)") => data.sediment_elevation[loc..., :]
         for (i, loc) in enumerate(grid_locations))...)
end

"""
    extract_sc(header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})

Extract Stratigraphic Column (SC) from the data. Returns a `DataFrame` with `time` and `sc<n>` columns where `<n>`
is in the range `1:length(grid_locations)`.
"""
function extract_sc(header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})
    n_facies = size(data.production)[1]
    DataFrame("time" => header.axes.t[1:end-1],
        ("sc$(i)_f$(f)" => stratigraphic_column(header, data, loc, f)
         for f in 1:n_facies, (i, loc) in enumerate(grid_locations))...)
end

"""
    extract_wd(header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})

Extract the water depth from the data. Returns a `DataFrame` with `time` and `wd<n>` columns where `<n>`
is in the range `1:length(grid_locations)`.
"""
function extract_wd(header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})
    na = [CartesianIndex()]
    wd = header.subsidence_rate .* header.axes.t[na, na, :] .- header.initial_topography[:, :, na] .- data.sediment_elevation
    DataFrame("time" => header.axes.t,
        ("wd$(i)" => wd[loc..., :] for (i, loc) in enumerate(grid_locations))...)
end
# ~/~ end

end
# ~/~ end
