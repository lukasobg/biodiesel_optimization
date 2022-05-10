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

# --------------------------- FINAL MODEL ------------------------------------
results = []

# set params
suppliers = [1000] 
capillarities = [10]

loop_times = []
t_main = @elapsed begin
    for sup in suppliers 
        for cap in capillarities
            Random.seed!(1234); # create same instance(s) each time
            for i in 1:50
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


res_file = string("../final_model/1000_10.csv")
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
