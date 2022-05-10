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

# --------------------------- Benchmarks ------------------------------------
t_total = @elapsed begin
    cap_factor = 2
    suppliers = 100

    file = "../benchmarks/$(suppliers)suppliers_$(cap_factor)capf.csv"

    bm_times = benchmark_model(suppliers, cap_factor)

    save_benchmark(bm_times, file)
end

println()
@printf "Total time taken: %.2f seconds \n" t_total
@printf "in minutes: %.2f \n" t_total/60
@printf "in hours: %.2f \n" t_total/3600
