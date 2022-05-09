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
using Printf
import Random

include("create_instances.jl")
include("run_models.jl")
include("benchmarks.jl")
include("save_results.jl")

Random.seed!(1234);

# --------------------------- Simple example problem ------------------------------------
#=
# Create the data for the example problem
toy_ins = create_toy_instance()

# Initialize models
m0_det = Model(Gurobi.Optimizer) # NMDT - deterministic
m0_nlm = Model(Gurobi.Optimizer) # nonlinear

# Create the models with the problem data from toy_ins
m_nlm = create_nonlinear_model(m0_nlm, toy_ins)

# Optimize the models with the problem data from toy_ins
println("\nNONLINEAR (TOY) MODEL OPTIMISATION\n")
m_nlm_opt = run_nonlinear_model(m_nlm) 

# Results
results = get_results(m_nlm_opt,"-","det")
=#
# --------------------------- Benchmarks ------------------------------------
#=
t_total = @elapsed begin
    cap_factor = 10
    suppliers = 250 

    file = "../benchmarks/$(suppliers)suppliers_$(cap_factor)capf.csv"
    println(file)

    bm_times = benchmark_model(suppliers, cap_factor)

    save_benchmark(bm_times, file)
end
minutes = t_total/60
println("Total time taken: $(t_total)")
println("in minutes: $(minutes)")
=#
# --------------------------- FINAL MODEL ------------------------------------
results = []

# set params
suppliers = [1000] 
capillarities = [6]

loop_times = []
t_main = @elapsed begin
    for sup in suppliers 
        for cap in capillarities
            Random.seed!(1234); # create same instance(s) each time
            for i in 1:6
                # create det instance
                t_loop = @elapsed begin
                    t_data = @elapsed begin
                        det_ins = create_deterministic_instance(sup,cap)
                    end

                    # DETERMINISTIC MODEL
                    # - init model
                    # - create model
                    # - optimise
                    # - save results
                    M0 =      Model(Gurobi.Optimizer)
                    M =       create_nonlinear_model(M0, det_ins)
                    M_opt =   run_nonlinear_model(M)
                    res_det = get_results(M_opt,sup,cap,"-")
                    push!(results, res_det)


                    # ROBUST MODELS WITH DIFF RISK LEVEL
                    # create rob instance
                    #=
                    data_entries = 100 # for SVC model in robust data generation
                    rob_ins = create_robust_instance(det_ins, data_entries)
                    risk_levels = [0.01,0.1,0.2]
                    for v in risk_levels
                        # - create svc ins
                        # - init model
                        # - create model
                        # - turn robust
                        # - optimise
                        # - save results
                        svc_ins = create_SVC_instance(rob_ins, v)
                        local M0 = Model(Gurobi.Optimizer)
                        local M = create_nonlinear_model(M0, det_ins)
                        M_rob,t = turn_model_robust(M, det_ins, rob_ins, svc_ins)
                        local M_opt = run_nonlinear_model(M_rob)
                        res_rob = get_results(M_opt,sup,cap,v)
                        push!(results, res_rob)
                        println("-------------------------------------------------------------------------- MODEL RISK $v SOLVED")
                    end
                    =#
                end
                push!(loop_times, t_loop)
                println("################### INSTANCE $(i) DONE, DATA GEN. $(round.(t_data/60; digits=2)), TOTAL $(round.(t_loop/60; digits=2)) MINUTES ###################")
                println()
            end
        end
    end
end


res_file = string("../final_model/test_run4.csv")
open(res_file, "a") do f
     write(f, "\n\n")
     write(f, "sup;cap;risk;PFAD;AF;CO;total biom;prod biod;demand;quality;prod cost;trans cost;process cost;hydro cost;total cost\n")
     writedlm(f, results, ';')
end

println()
println("Time taken in each instance:")
for i in 1:length(loop_times)
    print("$(i): ")
    @printf "%.2f sec, " loop_times[i]
    @printf "%.2f min" loop_times[i]/60
    println()
end

println()
@printf "Total time taken: %.2f seconds \n" t_main
@printf "in minutes: %.2f \n" t_main/60
@printf "in hours: %.2f \n" t_main/3600
