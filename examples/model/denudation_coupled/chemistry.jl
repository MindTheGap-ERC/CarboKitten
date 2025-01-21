using CSV, DataFrames, GLMakie

dis_adm = CSV.read("data/output/dissolution_adm.csv", DataFrame)
phys_adm = CSV.read("data/output/physical_adm.csv", DataFrame)
emp_adm = CSV.read("data/output/empirical_adm.csv", DataFrame)
nd_adm = CSV.read("data/output/nodenudation_adm.csv", DataFrame)
dis_adm_shallow = dis_adm[1:end,2]
phys_adm_shallow = phys_adm[1:end,2]
emp_adm_shallow = emp_adm[1:end,2]
nd_adm_shallow = nd_adm[1:end,2]
time = dis_adm[1:end,1]
d18O = -2 .* sin.(2π .* time / 0.2)

dis_preserved = diff(dis_adm_shallow)
preservation_ind = Array{Union{Nothing,Float64}}(undef,size(dis_preserved))
for i in 1:5000
    if  dis_preserved[i] >= 0.0001
        preservation_ind[i] = 1.0
    else 
        preservation_ind[i] = NaN 
    end
end
return preservation_ind
d18O_preserved = d18O[2:end] .* preservation_ind   
f = Figure()
ax = Axis(f[1, 1], xlabel = "Thickness (m)", ylabel = "δ¹⁸O (‰)")
scatter!(ax,dis_adm_shallow[2:end],d18O_preserved)


dis_adm_shallow[2:20:end]
d18O_sampled = d18O_preserved[1:20:end]
ax2 = Axis(f[1, 2], xlabel = "Time (Ma)", ylabel = "δ¹⁸O (‰)")

l1 = scatter!(ax2,time[2:end],d18O[2:end])
l2 = scatter!(ax2,time[2:end],d18O_preserved)
Legend(f[1,3],[l1,l2],["Original","Preserved"])
save("docs/src/_fig/chemistry.png",f)