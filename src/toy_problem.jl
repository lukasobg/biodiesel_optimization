using JuMP, Gurobi, Cbc, Ipopt
using LinearAlgebra, Statistics
using Distributions
using DelimitedFiles
using Printf
import Random

include("create_instances.jl")
include("run_models.jl")
include("benchmarks.jl")
include("save_results.jl")

Random.seed!(1234);

# --------------------------- Simple example problem ------------------------------------
# Create the data for the example problem
s = 12 # n of suppliers
toy_ins = create_toy_instance(s)

# Initialize models
m0_nlm = Model(Gurobi.Optimizer) # nonlinear

# Create the models with the problem data from toy_ins
m_nlm = create_nonlinear_model(m0_nlm, toy_ins)

# Optimize the models with the problem data from toy_ins
m_nlm_opt = run_nonlinear_model(m_nlm) 

# Results
println("---------------- NONLINEAR -----------------")
res_nlm = get_results(m_nlm_opt,s,0,"nlm")

res = []
push!(res, res_nlm)

res_file = string("../toy_problem/comparison.csv")
open(res_file, "a") do f
     write(f, "\n\n")
     write(f, "sup;cap;risk;PFAD;AF;CO;total biom;prod biod;demand;quality;prod cost;trans cost;process cost;hydro cost;total cost\n")
     writedlm(f, res, ';')
end

println("Biomasses before blending:")
println(value.(m_nlm_opt[:y]))

#=
# Initialize models
m0_lin = Model(Gurobi.Optimizer) # NMDT - deterministic

# Create the models with the problem data from toy_ins
m_lin = create_linear_model(m0_lin, toy_ins, 4)

# Optimize the models with the problem data from toy_ins
m_lin_opt = run_linear_model(m_lin) 

# Results
println("---------------- LINEAR NMDT -----------------")
res_lin = get_results(m_lin_opt,s,0,"lin-nmdt")
=#
