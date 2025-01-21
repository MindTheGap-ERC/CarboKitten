using CSV, DataFrames, GLMakie

dis_sac = CSV.read("data/output/dissolution_sac.csv", DataFrame)
phys_sac = CSV.read("data/output/physical_sac.csv", DataFrame)
emp_sac = CSV.read("data/output/empirical_sac.csv", DataFrame)
nd_sac = CSV.read("data/output/nodenudation_sac.csv", DataFrame)
time = dis_sac[1:end,1]
dis_sac_shallow = dis_sac[1:end,2]
phys_sac_shallow = phys_sac[1:end,2]
emp_sac_shallow = emp_sac[1:end,2]
nd_sac_shallow = nd_sac[1:end,2]
sac_total = hcat(time,dis_sac_shallow,phys_sac_shallow,emp_sac_shallow,nd_sac_shallow)
fig2 = Figure(resolution = (600, 800))
ax2 = Axis(fig2[1, 1],title = "Sediment accumulation Plots", xlabel = "Age (Myr)", ylabel = "Thickness (m)")
mode = ["Dissolution","Physical Erosion", "Empirical Denudation","No Denudation"]

for i in 1:4
    lines!(ax2, time, sac_total[:,i+1], label = mode[i])
end
fig2[1, 2] = Legend(fig2, ax2)
save("docs/src/_fig/sac_denudation.png", fig2)