# --------------------------------------- Setting up ------------------------------------
# Running for the first time:
#
# using Pkg
# Pkg.add("JuMP")
# Pkg.add("Cbc")
# Pkg.add("Ipopt")
# Pkg.add("Distributions")
# Pkg.add("DelimitedFiles")
#
# Pkg.add("Gurobi") # requires manual installation and a license
# Pkg.build("Gurobi")
#
# Gurobi requires separate installation and a license

# NOTE 1: check order of using, includes
using JuMP, Gurobi, Cbc, Ipopt
using LinearAlgebra, Statistics
using Distributions
using DelimitedFiles
import Random

include("create_instances.jl")
include("run_models.jl")
include("benchmarks.jl")
include("model_results.jl")

Random.seed!(1234);

# --------------------------- Simple example problem ------------------------------------
#= OUTDATED FUNCTION CALLS (PARAM CHANGES)

println("Small problem start")

# Create the data for the example problem
toy_ins = create_toy_instance()

# Initialize models
m0_det = Model(Gurobi.Optimizer) # NMDT - deterministic
m0_nlm = Model(Gurobi.Optimizer) # nonlinear

# Create the models with the problem data from toy_ins
m_det = create_linear_model(m0_det, toy_ins)
m_nlm = create_nonlinear_model(m0_nlm, toy_ins)

# Optimize the models with the problem data from toy_ins
println("\nNMDT DETERMINISTIC (TOY) MODEL OPTIMISATION\n")
m_det_opt = run_linear_model(m_det)
println("\nNONLINEAR (TOY) MODEL OPTIMISATION\n")
m_nlm_opt = run_nonlinear_model(m_nlm) 

# Results
println("\nNMDT DETERMINISTIC (TOY) MODEL RESULTS\n")
println("x_det: $(value.(m_det_opt[:x]))");
println("obj_det: $(objective_value(m_det_opt))");
println("solve_time: $(solve_time(m_det_opt))");

println("\nNONLINEAR (TOY) MODEL RESULTS\n")
println("x_det: $(value.(m_nlm_opt[:x]))");
println("obj_det: $(objective_value(m_nlm_opt))");
println("solve_time: $(solve_time(m_nlm_opt))");

=#
# --------------------------- Large example problem --------------------------------------
#= OUTDATED FUNCTION CALLS (PARAM CHANGES)

println("Large problem start")

# Create the data for the larger problem
suppliers = 20 
data_entries = 100 # for SVC model in robust data generation
det_ins, rob_ins = create_general_instance(suppliers, data_entries)

# Initialize models
#M0_det = Model(Gurobi.Optimizer) # NMDT - deterministic 
M0_rob = Model(Gurobi.Optimizer) # NMDT - robust 
#M0_nlm = Model(Gurobi.Optimizer) # nonlinear

# Create the models with the problem data from det_ins and rob_ins
#M_det = create_linear_model(M0_det, det_ins)

M_rob = create_linear_model(M0_rob, det_ins)
M_rob = turn_model_robust(M_rob, det_ins, rob_ins)

#M_nlm = create_nonlinear_model(M0_nlm, det_ins)

# Optimize the models with the problem data from det_ins and rob_ins
println("\nNMDT DETERMINISTIC MODEL OPTIMISATION\n")
M_det_opt = run_linear_model(M_det)
println("\nNMDT ROBUST MODEL OPTIMISATION\n")
M_rob_opt = run_linear_model(M_rob) 
println("\nNONLINEAR MODEL OPTIMISATION\n")
M_nlm_opt = run_nonlinear_model(M_nlm) 

# Results
println("\nNMDT DETERMINISTIC MODEL RESULTS\n")
println("x_det: $(value.(M_det_opt[:x]))");
println("obj_det: $(objective_value(M_det_opt))");
println("solve_time: $(solve_time(M_det_opt))");

println("\nNMDT ROBUST MODEL RESULTS\n")
println("x_det: $(value.(M_rob_opt[:x]))");
println("obj_det: $(objective_value(M_rob_opt))");
println("solve_time: $(solve_time(M_rob_opt))");
# @info value.(M_rob_opt[:slack]) == 0.0 ? "Model is feasible." : "Model is infeasible. Demand not met: $(value.(slack))"

println("\nNONLINEAR MODEL RESULTS\n")
println("x_det: $(value.(M_nlm_opt[:x]))");
println("obj_det: $(objective_value(M_nlm_opt))");
println("solve_time: $(solve_time(M_nlm_opt))");
=#

# --------------------------- Benchmarks ------------------------------------
t_total = @elapsed begin
    file = "../benchmarks/deterministic/______.csv"

    data_entries = 50 #100 # for SVC model in robust data generation
    suppliers = 100 

    bm_times = benchmark_model(data_entries, suppliers)

    save_benchmark(bm_times, file)
end
minutes = t_total/60
println("Total time taken: $(t_total)")
println("in minutes: $(minutes)")

# --------------------------- FINAL MODEL ------------------------------------
#=

# set params
suppliers = 50
data_entries = 50 # for SVC model in robust data generation
accuracy_p = 2

capillarity_factors = [2,3,5]
risk_levels = [0.01,0.1,0.2]

for c in capillarity_factors

    # Run robust models
    for v in risk_levels

        # create data instance
        det_ins, rob_ins = create_general_instance(suppliers, data_entries, c, v)

        # create model
        M0 = Model(Gurobi.Optimizer)
        M = create_linear_model(M0, det_ins, accuracy_p)
        M,t = turn_model_robust(M, det_ins, rob_ins)

        # optimise model
        M_opt = run_linear_model(M)

        # get results
        res_file = string("../final_model/model2/capf_",c,"_risk_",v,".csv")

        if result_count(M_opt) > 0
            get_results(M_opt,res_file)
        end

        println("-------------------------------------------------------------------------- MODEL CAP $c RISK $v SOLVED")

    end

    # Run deterministic model

    # create data instance
    det_ins, rob_ins = create_general_instance(suppliers, data_entries, c, 0.1)

    # create model
    M0 = Model(Gurobi.Optimizer)
    M = create_linear_model(M0, det_ins, accuracy_p)

    # optimise model
    M_opt = run_linear_model(M)

    # get results
    res_file = string("final_model/capf_",c,"_det.csv")
    
    if result_count(M_opt) > 0
        get_results(M_opt,res_file)
    end

end
=#
