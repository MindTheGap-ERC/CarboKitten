using CSV, DataFrames, GLMakie

dis_adm = CSV.read("data/output/dissolution_adm.csv", DataFrame)
phys_adm = CSV.read("data/output/physical_adm.csv", DataFrame)
emp_adm = CSV.read("data/output/empirical_adm.csv", DataFrame)
nd_adm = CSV.read("data/output/nodenudation_adm.csv", DataFrame)
time = dis_adm[1:end,1]
dis_adm_shallow = dis_adm[1:end,2]
phys_adm_shallow = phys_adm[1:end,2]
emp_adm_shallow = emp_adm[1:end,2]
nd_adm_shallow = nd_adm[1:end,2]

adm_total = hcat(time,dis_adm_shallow,phys_adm_shallow,emp_adm_shallow,nd_adm_shallow)