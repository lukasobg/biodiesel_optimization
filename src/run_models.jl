# Function creates a linear model m
#
# The data used in the variables, constraints and objective
# come from the problem instance parameter ins created in
# create_instances.jl
#
# Constraints #3 and #6 are replaced with NMDT
#
# Returns model m
function create_linear_model(m, ins)

    println("calling create_linear_model()")
    
    # ------------------ Get problem data from ins -------------------
    b1 = ins.b1; b = 1:b1 
    s1 = ins.s1; s = 1:s1
    k = ins.k
    p = ins.p; l = 1:p

    a_b = ins.a_b
    D_tot = ins.D_tot

    x_sb_min = ins.x_sb_min
    x_sb_max = ins.x_sb_max

    V_b_min = ins.V_b_min
    V_b_max = ins.V_b_max

    rho_sb = ins.rho_sb
    rho_min = ins.rho_min
    rho_max = ins.rho_max

    α = ins.α

    CB_b = ins.CB_b
    PC_b = ins.PC_b
    TC_sb = ins.TC_sb
    HC = ins.HC

    # -------------------- Set the problem variables -------------------
    # Basic variables (with domains) x, y, v, rho, q, z
    @variable(m, 0<= x[s,b]); #amount of biomass b from supplier s, tons 
    @variable(m, y[b] >= 0); #Amount of biomass b pre-treated for blending
    @variable(m, v[b] >= 0); #Property value of pooled biomass b
    @variable(m, rho_min <= rho <= rho_max); #Property value of blended biomass
    @variable(m, q >= 0); #Amount of bio-diesel before hydro pre-treatment
    @variable(m, z >= 0); #Amount of final bio-diesel
    # NMDT variables
    @variable(m, z1[k,l,b], Bin); # binary variable for v_b
    @variable(m, z2[k,l], Bin); # binary variable for rho 
    @variable(m, w1[s,b]); # replaces v_b * x_[s,b]
    @variable(m, Δw1[s,b]); # continuous domain for w1
    @variable(m, x1[k,l,s,b]); # disag. var. from lin. of x_[s,b]*z1_[k,l,s,b]
    @variable(m, Δv[b]>=0); # continuous domain of v_b in [0,10^{-p}]
    @variable(m, w2); # replaces p*q
    @variable(m, Δw2>=0); # continuous domain of w2
    @variable(m, x2[k,l] >=0); # disag. var. from lin. of q*z2_[k,l]
    @variable(m, Δρ>=0); # continuous domain of rho in [0,10^{-p}]

    # ------------------ Set the problem constraints --------------------
    # Basic constraints (2,4,5,8,9,10)
    @constraint(m, supply_limit[i=s, j=b], x_sb_min[i,j] <= x[i,j] <= x_sb_max[i,j]);  #2
    #constraint 3 handled with NMDT below
    @constraint(m, blendingrho[j=b], V_b_min[j] <= v[j] <= V_b_max[j]); #4 
    @constraint(m, pretreatment[j=b], y[j]==sum(a_b[j]*x[i,j] for i=1:s1)); #5
    #constraint 6 handled with NMDT below
    #constraint 7 (rho_limit) done at variable rho defenition
    @constraint(m, beforehydro, sum(y[j] for j=1:b1)==q); #8
    @constraint(m, hydro_conversion, α*q==z); #9
    @constraint(m, demand_fullfillment, z>=D_tot); #10

    # ---------------------- Set NMDT constraints -----------------------
    ########@constraint(m, eq11[i=s, j=b], w1[i,j]==v[j]*x[i,j]);
    @constraint(m, eq12[i=s, j=b], w1[i,j]==sum(sum( 10.0^(-l)*k*x1[k,l,i,j] for l=1:p) for k=0:9)+Δw1[i,j]);
    @constraint(m, eq13[j=b], v[j]==sum(sum( 10.0^(-l)*k*z1[k,l,j] for l=1:p) for k=0:9)+Δv[j]);

    @constraint(m, eq14[i=s, j=b, l=1:p], x[i,j]==sum(x1[k,l,i,j] for k=0:9));

    @constraint(m, eq151[k=k,l=1:p, i=s, j=b], x1[k,l,i,j] <= x_sb_max[i,j]*z1[k,l,j]) ;
    @constraint(m, eq152[k=k,l=1:p, i=s, j=b], x_sb_min[i,j]*z1[k,l,j] <= x1[k,l,i,j]) ; 

    @constraint(m, eq16[l=1:p, j=b], sum( z1[k,l,j] for k=0:9)==1);

    @constraint(m, eq171[i=s, j=b], Δw1[i,j]<= x_sb_max[i,j]*Δv[j]);
    @constraint(m, eq172[i=s, j=b], x_sb_min[i,j]*Δv[j] <= Δw1[i,j]);

    @constraint(m, eq18[i=s, j=b], Δw1[i,j]<= (x[i,j] - x_sb_min[i,j])*10.0^(-p)+x_sb_min[i,j]*Δv[j]);
    @constraint(m, eq19[i=s, j=b], Δw1[i,j]>= (x[i,j] - x_sb_max[i,j])*10.0^(-p)+x_sb_max[i,j]*Δv[j]);
    @constraint(m, eq20[j=b], Δv[j]<= 10.0^(-p));
    #########@constraint(m, eq21, w2== rho*q); #not in the model
    @constraint(m, eq22, w2== sum(sum( 10.0^(-l)*k*x2[k,l] for l=1:p) for k=0:9) + Δw2);
    @constraint(m, eq23, rho== sum(sum( 10.0^(-l)*k*z2[k,l] for l=1:p) for k=0:9) + Δρ);

    @constraint(m, eq24[l=1:p], q== sum(x2[k,l] for k=0:9)); 

    @constraint(m, eq25[k=0:9, l=1:p], x2[k,l] <= sum(sum(x_sb_max[i,j] for j=1:b1) for i=1:s1)* z2[k,l]);
    @constraint(m, eq26[l=1:p], sum( z2[k,l] for k=0:9)==1);
    @constraint(m, eq27, Δw2<= sum(sum(x_sb_max[i,j] for j=1:b1) for i=1:s1)* Δρ);
    @constraint(m, eq28, Δw2<= (q-0)*10.0^(-p)+0*Δρ);
    @constraint(m, eq29, Δw2>= (q-sum(sum(x_sb_max[i,j] for j=1:b1) for i=1:s1))*10.0^(-p)+sum(sum(x_sb_max[i,j] for j=1:b1) for i=1:s1)*Δρ);
    @constraint(m, eq30, Δρ<= 10.0^(-p));
    @constraint(m, eq31[j=1:b1], sum(rho_sb[i]*x[i,j] for i=1:s1)==sum(w1[i,j] for i=1:s1));
    @constraint(m, eq32, w2==sum(sum(a_b[j]*w1[i,j]  for i=1:s1) for j=1:b1));

    
    # ------------------------ Set the objective --------------------------
    @objective(m, Min, sum(sum((CB_b[j]+TC_sb[i]+PC_b[j])*x[i,j] for j=1:b1) for i=1:s1)+HC*q);

    return m
end

# Function optimises the linear model m
#
# Returns optimsed model m
function run_linear_model(m)

    println("calling run_linear_model()")

    # Optimize and return the optimized model m
    println()
    println("------------- Output from !optimize(m) call -------------")
    println()

    optimize!(m);

    println()
    println("---------- End of output from !optimize(m) call ---------")
    println()
    
    return m
end

# Function creates nonlinear model m
#
# The data used in the variables, constraints and objective
# come from the problem instance parameter ins
#
# Constraints #3 and #6 are INCLUDED -> nonlinear optimisation
#
# Used for validation of results
#
# Returns model m
function create_nonlinear_model(m, ins)

    println("calling create_nonlinear_model()")
    
    # ------------------ Get problem data from ins -------------------
    b1 = ins.b1; b = 1:b1 
    s1 = ins.s1; s = 1:s1
    #k = ins.k # NMDT
    #p = ins.p; l = 1:p # NMDT

    a_b = ins.a_b
    D_tot = ins.D_tot

    x_sb_min = ins.x_sb_min
    x_sb_max = ins.x_sb_max

    V_b_min = ins.V_b_min
    V_b_max = ins.V_b_max

    rho_sb = ins.rho_sb
    rho_min = ins.rho_min
    rho_max = ins.rho_max

    α = ins.α

    CB_b = ins.CB_b
    PC_b = ins.PC_b
    TC_sb = ins.TC_sb
    HC = ins.HC

    # -------------------- Set the problem variables -------------------
    # Variable domains for x, y, v, rho, q, z
    @variable(m, 0<= x[s,b]); #amount of biomass b from supplier s, tons 
    @variable(m, y[b] >= 0); #Amount of biomass b pre-treated for blending
    @variable(m, v[b] >= 0); #Property value of pooled biomass b
    @variable(m, rho_min <= rho <= rho_max); #Property value of blended biomass
    @variable(m, q >= 0); #Amount of bio-diesel before hydro pre-treatment
    @variable(m, z >= 0); #Amount of final bio-diesel
    
    # ------------------ Set the problem constraints --------------------
    # Create basic constraints (2,3,4,5,6,8,9,10)
    @constraint(m, supply_limit[i=s, j=b], x_sb_min[i,j] <= x[i,j] <= x_sb_max[i,j]);  #2
    @constraint(m, blending1[j=b], sum(x[i,j]*rho_sb[i] for i=1:s1)==sum(x[i,j]*v[j] for i=1:s1)) #3 
    @constraint(m, blendingrho[j=b], V_b_min[j] <= v[j] <= V_b_max[j]); #4 
    @constraint(m, pretreatment[j=b], y[j]==sum(a_b[j]*x[i,j] for i=1:s1)); #5
    @constraint(m, blending2, sum( v[j]*y[j] for j=1:b1)==rho*q) #6
    #constraint 7 (rho_limit) done at variable rho defenition
    @constraint(m, beforehydro, sum(y[j] for j=1:b1)==q); #8
    @constraint(m, hydro_conversion, α*q==z); #9
    @constraint(m, demand_fullfillment, z>=D_tot); #10
     
    # ------------------------ Set the objective --------------------------
    @objective(m, Min, sum(sum((CB_b[j]+TC_sb[i]+PC_b[j])*x[i,j] for j=1:b1) for i=1:s1)+HC*q);

    return m
end

# Function optimises the nonlinear model m
#
# Used for validation
#
# Returns optimsed model m
function run_nonlinear_model(m) 

    println("calling run_nonlinear_model()")

    # ------------ Configure optimizer for non-linear model ---------------
    set_optimizer_attributes(m,
        "NonConvex" => 2,
        "NumericFocus" => 3,
        "OptimalityTol" => 1e-9,
        "MIPGap" => 1e-9,
        "FeasibilityTol" => 1e-9,
        "IntFeasTol" => 1e-9,
        "Presolve" => 1,
        "ScaleFlag" => 3
    )

    # Optimize and return the optimized model m
    println()
    println("------------- Output from !optimize(m) call -------------")
    println()

    optimize!(m);

    println()
    println("---------- End of output from !optimize(m) call ---------")
    println()

    return m
end

# Function makes the model m robust
# 
# Adds robust variables and constraints
# Updates the objective

# The data used in the variables, constraints and objective
# come from the problem instance parameters det_ins and rob_ins 
# created in create_instances.jl
#
# Returns robust version of model m
function turn_model_robust(m, det_ins, rob_ins)

    println("calling turn_model_robust()")

    # ------------------ Get problem data from det_ins -------------------
    b1 = det_ins.b1; b = 1:b1
    s1 = det_ins.s1; s = 1:s1
    #k = det_ins.k # NMDT
    #p = ins.p; l = 1:p # NMDT

    a_b = det_ins.a_b 
    D_tot = det_ins.D_tot

    #x_sb_min = det_ins.x_sb_min
    #x_sb_max = det_ins.x_sb_max

    #V_b_min = det_ins.V_b_min
    #V_b_max = det_ins.V_b_max

    #rho_sb = det_ins.rho_sb
    #rho_min = det_ins.rho_min
    #rho_max = det_ins.rho_max

    α = det_ins.α

    # Needed only if objective is updated
    CB_b = det_ins.CB_b
    TC_sb = det_ins.TC_sb
    PC_b = det_ins.PC_b
    HC = det_ins.HC

    # ------------------ Get problem data from rob_ins -------------------
    B_sb = rob_ins.B_sb

    D = rob_ins.D
    Q = rob_ins.Q
    θ = rob_ins.θ
    α_val = rob_ins.α_val

    set_SV = rob_ins.set_SV
    set_SVB = rob_ins.set_SVB

    # Variables from deterministic used in constraints here
    x = m[:x]
    q = m[:q] # only needed if objective updated here

    # Robust variables
    @variable(m, λ[j in s, i in set_SV] >= 0);
    @variable(m, μ[j in s, i in set_SV] >= 0);
    @variable(m, η >= 0);
    @variable(m, slack >= 0);

    # Robust constraints

    println()
    println("Time taken in constraint 41)")
    @time begin

    # 41 first alpha is alpha from det model
    # (processing yield of biomass b)
    @constraint(m, [ii in set_SVB],
        α * (
            sum(a_b[ib] *
                sum(
                    sum(D[s2,iSV] * B_sb[s2,ib] * # UNSURE: B_sb needed here or not?
                        sum(
                            (λ[s1,iSV]-μ[s1,iSV]) * Q[s1,s2] 
                        for s1 in s)
                    for s2 in s)
                for iSV in set_SV) 
            for ib in b)) - 
        η*θ[ii] >= D_tot * 100 - slack # UNSURE: use slack or not
    );

    end
    println()

    # 42 
    @constraint(m, [j in s, i in set_SV],
        λ[j,i] + μ[j,i] - η * α_val[i] == 0
    )

    # 43
    @constraint(m, [k in b, j in s],
        sum(sum(Q[j,j1] * (λ[j1,i] - μ[j1,i]) for j1 in s) for i in set_SV) == x[j,k]
    )

    # UNSURE: Update objective to have slack or not 
    @objective(m, Min,
        sum(
             sum(
                 (CB_b[j]+TC_sb[i]+PC_b[j])*x[i,j] 
             for j=1:b1)
         for i=1:s1) + 
         HC*q + 1e3*slack
    );

    return m
end
