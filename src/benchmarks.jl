# If running deterministic:
#  - comment out turn_model_robust() calls in # Build section
#  - return times WITHOUT t_cons times
#
# If running robust:
#  - comment in turn_model_robust() calls in # Build section
#  - return times WITH t_cons times
function benchmark_model(suppliers, cap_factor)

    times = []

    for i in 1:11
        println("-------------------------------------------------------------------------- INSTANCE $i START")
        # Create the data for the larger problem

        det_ins = create_deterministic_instance(suppliers, cap_factor)

        M0_nlm = Model(Gurobi.Optimizer)
        M0_nmdt2 = Model(Gurobi.Optimizer)
        M0_nmdt3 = Model(Gurobi.Optimizer)
        M0_nmdt4 = Model(Gurobi.Optimizer)

        ###### Build and time models

        t_build1 = @elapsed begin
            M_nlm = create_nonlinear_model(M0_nlm, det_ins)
        end
        println("-------------------------------------------------------------------------- INSTANCE $i MODEL 1 BUILT")

        t_build2 = @elapsed begin
            M_nmdt2 = create_linear_model(M0_nmdt2, det_ins, 2)
        end
        println("-------------------------------------------------------------------------- INSTANCE $i MODEL 2 BUILT")

        t_build3 = @elapsed begin
            M_nmdt3 = create_linear_model(M0_nmdt3, det_ins, 3)
        end
        println("-------------------------------------------------------------------------- INSTANCE $i MODEL 3 BUILT")

        t_build4 = @elapsed begin
            M_nmdt4 = create_linear_model(M0_nmdt4, det_ins, 4)
        end
        println("-------------------------------------------------------------------------- INSTANCE $i MODEL 4 BUILT")

        ###### Solve and time models
        t_solve1 = @elapsed begin
            M_nlm_opt = run_nonlinear_model(M_nlm)
        end
        println("-------------------------------------------------------------------------- INSTANCE $i MODEL 1 SOLVED")

        t_solve2 = @elapsed begin
            M_nmdt2_opt = run_linear_model(M_nmdt2)
        end
        println("-------------------------------------------------------------------------- INSTANCE $i MODEL 2 SOLVED")

        t_solve3 = @elapsed begin
            M_nmdt3_opt = run_linear_model(M_nmdt3)
        end
        println("-------------------------------------------------------------------------- INSTANCE $i MODEL 3 SOLVED")

        t_solve4 = @elapsed begin
            M_nmdt4_opt = run_linear_model(M_nmdt4)
        end
        println("-------------------------------------------------------------------------- INSTANCE $i MODEL 4 SOLVED")

        ###### Return times

        # deterministic
        push!(times, (t_build1,t_solve1,t_build2,t_solve2,t_build3,t_solve3,t_build4,t_solve4))

        #robust
        #push!(times, (t_build1,t_cons1,t_solve1,t_build2,t_cons2,t_solve2,t_build3,t_cons3,t_solve3,t_build4,t_cons4,t_solve4))
    end

    return times
end

# If running deterministic:
#  - write header WITHOUT cons times
#
# If running robust:
#  - write header WITH cons times
function save_benchmark(times, file)

    open(file; write=true) do f
         # deterministic
         write(f, "nlm build;nlm solve;nmdt2 build;nmdt2 solve;nmdt3 build;nmdt3 solve;nmdt4 build;nmdt4 solve\n")

         #robust
         #write(f, "nlm build;nlm cons;nlm solve;nmdt2 build;nmdt2 cons;nmdt2 solve;nmdt3 build;nmdt3 cons;nmdt3 solve;nmdt4 build;nmdt4 cons;nmdt4 solve\n")

         writedlm(f, times, ';')
    end

end
