# --------------------------------------- Setting up ------------------------------------
# Running for the first time:
#
# using Pkg
# Pkg.add("JuMP")
# Pkg.add("Cbc")
# Pkg.add("Ipopt")
# Pkg.add("Distributions")
#
# Pkg.add("Gurobi") # requires manual installation and a license
# Pkg.build("Gurobi")
#
# Gurobi requires separate installation and a license

# NOTE 1: check order of using, includes
using JuMP, Gurobi, Cbc, Ipopt
using LinearAlgebra, Statistics
using Distributions
import Random

include("create_instances.jl")
include("run_models.jl")

Random.seed!(1234);

# --------------------------- Simple example problem ------------------------------------

#=
println("Small problem start")

# Create the data for the example problem
ins = create_toy_instance()

# Model for deterministic 
m = Model(Gurobi.Optimizer)

# Optimize the model with the problem data from ins
det_model = run_linear_model(m, ins) # Deterministic

# Print results
#x_det = value.(det_model[:x1]);
#obj_det = objective_value(det_model);
#solve_time(det_model);

# Validate results with nonlinear deterministic model
# NOTE 2: should a new identical ins be created here?
#ins_nlm = create_toy_instance()
nlm = Model(Gurobi.Optimizer)
nlm_model = run_nonlinear_model(nlm, ins)
=#

# --------------------------- Larger example problem --------------------------------------

println("Large problem start")

# Create the data for the larger problem
det_ins, rob_ins = create_general_instance() 

#=
# Models for deterministic and robust
L_m = Model(Gurobi.Optimizer)
L_rm = Model(Gurobi.Optimizer)

# Optimize the models with the problem data from ins
L_det_model = run_linear_model(L_m, det_ins) # Deterministic
L_rob_model = run_robust_model(L_rm, det_ins, rob_ins) # Robust

# Print results

# Validate results with nonlinear deterministic model
# NOTE 2: should a new identical ins be created here?
#det_ins, rob_ins = create_general_instance() 
L_nlm = Model(Gurobi.Optimizer)
L_nlm_model = run_nonlinear_model(det_ins)
=#
