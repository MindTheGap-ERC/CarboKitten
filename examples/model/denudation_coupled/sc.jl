using CSV, DataFrames, GLMakie

data = CSV.read("data/output/dissolution_sc.csv", DataFrame)

time = data[:,1]
litho1 = data[:,2]
litho2 = data[:,3]
litho3 = data[:,4]
sum_facies = litho1 .+ litho2 .+ litho3
sum_facies = vcat(0,sum_facies)
fig =Figure()
ax = Axis(fig[1, 1])

x_left = 0.2  
x_right = 0.35

for i in 1:length(time)-1
    boundary1 = zeros(Float64,length(time))
    boundary2 = zeros(Float64,length(time))
    boundary3 = zeros(Float64,length(time))

    boundary1[i+1] = sum_facies[i] .+ boundary1[i] # bot of Lithology 1
    boundary2[i+1] = sum_facies[i] .+ boundary1[i+1]  # bot of Lithology 2
    boundary3[i+1] = sum_facies[i] .+ boundary2[i+1] # bot of Lithology 3
    hspan!(ax,boundary1[i], boundary2[i+1],color = :green)
    hspan!(ax, boundary2[i], boundary3[i+1], color = :red)
    hspan!(ax, boundary3[i], boundary3[i+1],color = :blue)
end

fig
save("docs/src/_fig/dissolution_sc.png", fig)