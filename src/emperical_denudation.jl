## this code used the emperical equations with form of D = f(precipitation, slope) to estimate the denudation rates on exposed carbonate platform.
module emperical_denudation

using CarboKitten.CaProd
export calculate_D
# calculate planar slopes based on [ARCGIS apporach](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-slope-works.htm)


function calculate_slope(elevation::Matrix{Float64}, cellsize::Float64) 
    nrows, ncols = size(elevation)
    slope = similar(elevation)

    padded_elevation = zeros(Float64, nrows + 2, ncols + 2)
    padded_elevation[2:end-1, 2:end-1] = elevation

    for i in 2:nrows-1
        for j in 2:ncols-1
            dzdx = ((padded_elevation[i - 1, j + 1] + 2*padded_elevation[i, j + 1] + padded_elevation[i + 1, j + 1]) - (padded_elevation[i - 1, j - 1] + 2*padded_elevation[i, j - 1] + padded_elevation[i + 1, j - 1])) ./ (8 * cellsize)
            dzdy = ((padded_elevation[i + 1, j - 1] + 2*padded_elevation[i + 1, j] + padded_elevation[i + 1, j + 1]) - (padded_elevation[i - 1, j - 1] + 2*padded_elevation[i - 1, j] + padded_elevation[i - 1, j + 1]))/ (8 * cellsize)
            slope[i, j] = atan.(sqrt.(dzdx.^2 + dzdy.^2))  * (180 / π)
        end
    end

    return slope
end


function calculate_D(precip::Float64, elevation::Matrix{Float64}, cellsize::Float64)
    slope = calculate_slope(elevation,cellsize)
    nrows, ncols = size(elevation)
    # function format
    for i in 1:nrows
        for j in 1:ncols
    D[i,j] = (9.1363 ./ (1 .+ exp.(-0.008519.*(precip .- 580.51)))) .* (9.0156 ./ (1 .+ exp.(-0.1245.*(slope[i,j] .- 4.91086)))) # using logistic function
        end
    end
    return D./1000 #m/kyr
end

end