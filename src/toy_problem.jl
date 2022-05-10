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
m0_lin = Model(Gurobi.Optimizer) # NMDT - deterministic
m0_nlm = Model(Gurobi.Optimizer) # nonlinear

# Create the models with the problem data from toy_ins
m_lin = create_linear_model(m0_lin, toy_ins, 3)
m_nlm = create_nonlinear_model(m0_nlm, toy_ins)

# Optimize the models with the problem data from toy_ins
m_lin_opt = run_linear_model(m_lin) 
m_nlm_opt = run_nonlinear_model(m_nlm) 

# Results
println("---------------- LINEAR NMDT -----------------")
res_lin = get_results(m_lin_opt,s,0,"lin-nmdt")
println("---------------- NONLINEAR -----------------")
res_nlm = get_results(m_nlm_opt,s,0,"nlm")
