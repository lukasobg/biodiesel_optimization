using JuMP, Gurobi, Cbc, Ipopt
using LinearAlgebra, Statistics
using Distributions
import Random

Random.seed!(1234);

m = Model(Gurobi.Optimizer)
#set_optimizer_attributes(m, "NonConvex" => 2,"Presolve"=>0)

#for sampling supplys
sdPFAD=11.5; #palm oil analytics database
avgPFAD=76;
sdAF=100; #European comission trade database , EU comission report of emissions of animal fats used for biodiesel 
avgAF=1500;
sdUCO=100; ##European comission trade database and UK biofuel statistics
avgUCO=500;

function create_supply(lambda, s, D_tot, sdPFAD, sdAF, sdUCO, avgPFAD, avgAF, avgUCO)

    PFAD=rand(DiscreteUniform(10, 15)) 
    WCO=rand(DiscreteUniform(20, 25))
    AF=s-PFAD-WCO
    Ds_test=zeros(Float64, s)
    PC_test=zeros(Float64, s)

    ######################### creation of B_sb
    B_sb_test=zeros(Int64, s, 3)
    for i in 1:PFAD
        B_sb_test[i,:]=[0 0 1] #PFAD
    end

    for k in PFAD+1:AF+PFAD
        B_sb_test[k,:]=[0 1 0] #AF
    end

    for j in AF+PFAD+1:s
        B_sb_test[j,:]=[1 0 0] #UCO
    end
    #########################
    
    ########################creation of distances
    for i in 1:PFAD
        dist_sample=rand(Normal(1100, 150 )) #*10 THESE ARE ONE WAY, move numbers out of the function?
        if (dist_sample<100)
            dist_sample=100
        end
        Ds_test[i]=dist_sample
    end

    for k in PFAD+1:PFAD+AF
        dist_sample=rand(Normal(330, 100)) #*10
        if (dist_sample<10)
            dist_sample=10
        end
        Ds_test[k]=dist_sample
    end

    for j in PFAD+AF+1:s
        dist_sample=rand(Normal(60, 30)) #*10
        if (dist_sample<2)
            dist_sample=2
        end
        Ds_test[j]=dist_sample
    end
    #########################
    

    #l=2.6 #to easily change sampling volumes
    ######################### creation of PC_s
    while( sum(PC_test) < D_tot )
       #PC_test=zeros(Float64, s)
       #Ds_test=zeros(Float64, s)
    for i in 1:PFAD
        sample=rand(Normal(avgPFAD, sdPFAD)) 
        if (sample<300)
            sample=300
        end
        PC_test[i]=sample
    end

    for k in PFAD+1:PFAD+AF
        sample=rand(Normal(avgAF, sdAF)) 
        if (sample<100)
            sample=100
        end
        PC_test[k]=sample
    end
    
    for j in PFAD+AF+1:s 
            sample=rand(Normal(avgUCO, sdUCO)) 
            if (sample<50)
                sample=50
            end
            PC_test[j]=sample
    end
    end

    return round.(PC_test, digits=2), round.(B_sb_test, digits=2), round.(Ds_test, digits=2)
end

#sets
b1=3; #number of biomass  
s1=12;  #number of suppliers 
b=1:b1; #types of biomass
s=1:s1; #suppliers
k=0:9;
p=3; #cardinal of l:1,2 or 3
l=1:p;

#a_b = [0.95 ;0.9; 0.99]; #oil percentages of biomass WCO, Animal fat, PFADs, old ones
a_b=[0.94; 0.95 ; 0.98];
p_low=0.05 #for big suppliers, for capillarity expansion
p_high=0.2 #for small suppliers
#D_tot = 3000; #*1000 #total biodiesel demand in tons
D_tot=1500;
lambda = 1.5;


PC_s, B_sb, D_s = create_supply(lambda, s1, D_tot, sdPFAD, sdAF, sdUCO, avgPFAD, avgAF, avgUCO)

B_sb=[0 0 1;0 0 1;0 0 1;0 0 1;0 1 0;0 1 0;0 1 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0]; #toy model test values

capillarity_factor=2; #give as input, test different levels

#CAHNGE: splitting randomly not to even parts
function create_capillarity_expansion(capillarity_factor, PC_s, B_sb,D_s, s1, r_low=50, r_high=100)
    UCO_suppliers=sum(B_sb[:,1])
    UCO_suppliers=Int(UCO_suppliers)
    #how many more UCO will be, all suppliers divided to #capillarity_factor many suppliers
    UCO_new=capillarity_factor*UCO_suppliers
    newB_sb=B_sb[1:(s1-UCO_suppliers),:]

    matrix=zeros(Int8, UCO_new , 3)
    matrix[:,1].=1
    newB_sb=vcat(newB_sb,matrix)
    #now go through PC_s of UCO suppliers and explode them to little suppliers
    copy=PC_s[s1-UCO_suppliers+1:s1]
    newPC_s=PC_s[1:s1-UCO_suppliers]
    for j in 1:length(copy)
        for k in 1:capillarity_factor
            value=copy[j]/capillarity_factor # randomdy split the capacity, five points
            append!(newPC_s, value)
        end
    end

    #Distanceses
    DsCopy=D_s[s1-UCO_suppliers+1:s1] #only uco distances
    newD_s=D_s[1:s1-UCO_suppliers]    #no uco distances
    
    for a in 1:length(DsCopy)
        for b in 1:capillarity_factor
            radius=rand(DiscreteUniform(r_low, r_high))
            val=DsCopy[a]+radius
            append!(newD_s, val)
        end
    end

    new_s1=s1-UCO_suppliers+UCO_new
    new_s=1:new_s1
    return newB_sb, newPC_s, newD_s, new_s1, new_s
end

B_sb, PC_s, D_s, s1, s= create_capillarity_expansion(capillarity_factor, PC_s, B_sb,D_s, s1)


function create_qualities(p_high, p_low, PC_s,B_sb,a_b,s1)
    interval=(p_high-p_low)/(s1-1)
    q_s=zeros(Float64,s1) #quality dimension
    for i in s
        j=1
        while PC_s[i] != (sort(PC_s))[j]
            j=j+1
            #println(j)
        end
        #println(j)
        index= j
        q_s[i]= rand(Normal((B_sb*a_b)[i], (p_high - (index-1)*interval)*(B_sb*a_b)[i] )) 
        
        if q_s[i]>1
            q_s[i]=1
        end
    end
    return q_s
end


rho_sb= create_qualities(p_high, p_low, PC_s,B_sb,a_b,s1) #rho_sb The property value of biomass bought from supplier s
rho_sb=[0.8 ;0.9; 0.78;0.64;0.68;0.85;0.7;0.58;0.7;0.65;0.55;0.6]; #toy model test values


x_sb_min=[20;30;12;10;10;20;10;15;10;10;10;10]; #toy model test values
x_sb_min=x_sb_min.*B_sb 

x_sb_max=[200;300;120;100;100;200;400;150;100;180;60;150]; #toy model test values
x_sb_max=x_sb_max.*B_sb
V_b_min=[0.6 0.6 0.6]; #Minimum property value of biomass b
V_b_max=[0.85 0.95 0.9]; #Maximum property value of biomass b

rho_min=0.6; #Minimum property value of bio-diesel
rho_max=0.9; #Maximum property value of bio-diesel
#Beta_b = [0.86;0.78;0.94]; #blending  yeilds, not in the new model
α=0.94; #hydroconversion yield?

#CB_b = [250, 400, 450]; #cost of biomass b, WCO , AF, PFAD , old ones
CB_b = [90, 120, 100];
#PC_b = [680 696 700]; #processing cost for biomass b, old ones   
PC_b = [50 60 40];
#T_b = 123; #transportation cost for biomass b  tons
#TC_sb=D_s*T_b #distance*transportation costs
TC_sb=[100;90;80;50;100;120;100;120;120;100;80;50]; #toy model test values
HC = 100; #hydrotreatment cost 

@variable(m, 0<= x[s,b]); #amount of biomass b from supplier s, tons 
@variable(m, y[b] >= 0); #Amount of biomass b pre-treated for blending
@variable(m, v[b] >= 0); #Property value of pooled biomass b
@variable(m, rho_min <= rho <= rho_max); #Property value of blended biomass
@variable(m, q >= 0); #Amount of bio-diesel before hydro pre-treatment
@variable(m, z >= 0); #Amount of final bio-diesel

##############################new ones###############################
@variable(m, z1[k,l,b], Bin); 
@variable(m, z2[k,l], Bin);
@variable(m, w1[s,b]); 
@variable(m, Δw1[s,b]); 
@variable(m, x1[k,l,s,b]); 
@variable(m, Δv[b]>=0);
@variable(m, w2);
@variable(m, Δw2>=0);
@variable(m, x2[k,l] >=0);
@variable(m, Δρ>=0);
##############################new ones###############################


function create_constraints()

    @constraint(m, supply_limit[i=s, j=b], x_sb_min[i,j] <= x[i,j] <= x_sb_max[i,j]);  #2

    #######@constraint(m, blending1[j=b], sum(x[i,j]*rho_sb[i] for i=1:s1)==sum(x[i,j]*v[j] for i=1:s1)) #3 not in model
    @constraint(m, blendingrho[j=b], V_b_min[j] <= v[j] <= V_b_max[j]); #4 
    @constraint(m, pretreatment[j=b], y[j]==sum(a_b[j]*x[i,j] for i=1:s1)); #5
    ########@constraint(m, blending2, sum( v[j]*y[j] for j=1:b1)==rho*q) #6 not in model

    #constraint 7 (rho_limit) done at variable rho defenition
    @constraint(m, beforehydro, sum(y[j] for j=1:b1)==q); #8
    @constraint(m, hydro_conversion, α*q==z); #9
    @constraint(m, demand_fullfillment, z>=D_tot); #10
    ##############################new ones###############################
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
end

create_constraints()



@objective(m, Min, sum(sum((CB_b[j]+TC_sb[i]+PC_b[j])*x[i,j] for j=1:b1) for i=1:s1)+HC*q);
println(m)

#=
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
=# 
 
optimize!(m)
x_det = value.(x)
obj_det = objective_value(m)

#######################END OF DETRMINISTIC######################

#data creation for robust model#
#gamma_s [0,1]
#gamma matrix Sx100 so past data for all suppliers
ent = 500;
D = zeros(Float64, length(s), ent); #saving gamma data here

for i in s
    if B_sb[i,1]==1
        Gamma_b=round.(rand(Normal(75, 5), ent), digits=2) #WCO
        D[i,:]=Gamma_b

    elseif B_sb[i,2]==1
        Gamma_b=round.(rand(Normal(95, 5), ent), digits=2)   #AF
        D[i,:]=Gamma_b

    else B_sb[i,3]==1
        Gamma_b=round.(rand(Normal(85, 5), ent), digits=2)   #PFAD
        D[i,:]=Gamma_b
    end
end


for i in s
    for j in 1:ent
        D[i,j] = (D[i,j]>100.00) ? 100.00 : D[i,j]
        D[i,j] = (D[i,j]<0.00) ? 0.00 : D[i,j]
    end
end

#D is now gamma_s s x entities
Σ = cov(D, corrected=true, dims=2);
Q = round.(Σ^(-0.5), digits = 6);


# Step 2: form the Kernel matrix
#J=1:1000
uU = maximum(D, dims=2); #dims=2 since data is per column
uL = minimum(D, dims=2);
N = s1; #N is the same as the number of suppliers
J = 1:N ;
K(u,v) = sum(uU[j] - uL[j] - abs(u[j] - v[j]) for j in J);

# The kernel trick K sweeps the data to be calculated.
M = ent; #data entities
Km = zeros(M, M);
for i = 1:M
    for i1 = 1:M
        Km[i,i1] = K(D[:,i], D[:,i1])
    end
end


# Step 3: derive SVC model
ν = 0.05;   # regularisation parameter
SVC_Model = Model(Ipopt.Optimizer);
set_optimizer_attributes(SVC_Model, "mehrotra_algorithm" => "yes");
@variable(SVC_Model, α[i in 1:M] >= 0);
@objective(SVC_Model, Min,
    sum(α[i] * α[j] * Km[i,j] for i in 1:M for j in 1:M)
    - sum(α[i] * Km[i,i] for i in 1:M)
);
@constraint(SVC_Model, [i in 1:M], α[i] * (M * ν) <= 1);
@constraint(SVC_Model, sum(α[i] for i in 1:M) == 1);

optimize!(SVC_Model)
α_val = round.(value.(α), digits = 6);

for i in 1:length(α_val)
    α_val[i] < 1e-6 ? α_val[i] = 0.0 : α_val[i]
end


# Filtering the support vectors and the boundary support vectors
ϵ = 1e-12 # tolerance for zero: helps finding boundary points correctly.
SV = [α_val[i] > 0 + ϵ ? 1 : 0 for i in 1:M]
SVB = [α_val[i] > 0 + ϵ  && α_val[i] < 1/(M * ν) - ϵ ? 1 : 0 for i in 1:M]

# Forming the sets with (boundary) support vectors
set_SV = []
set_SVB = []
for i in 1:M
    if SV[i] == 1
        push!(set_SV, i)
    end
    if SVB[i] == 1
        push!(set_SVB, i)
    end
end

θ = Dict();

for i1 in set_SVB
    θ[i1] = sum(α_val[i] * norm(Q * (D[:,i1] - D[:,i]), 1) for i in set_SV)
end

RobModel = Model(Gurobi.Optimizer)
@variable(RobModel, x[s] >= 0); #amount of biomass b from supplier s
@variable(RobModel, w >= 0);
@variable(RobModel, 0.78 <= rho <= 0.94);
@variable(RobModel, y[b] >= 0);
@variable(RobModel, z >= 0) ;

@variable(RobModel, λ[j in s, i in set_SV] >= 0);
@variable(RobModel, μ[j in s, i in set_SV] >= 0);
@variable(RobModel, η >= 0);
@variable(RobModel, slack >= 0);

@objective(RobModel, Min, sum(sum((C_b[i]+T_b[i]+D_s[j]+P_b[i])B_sb[j,i]*x[j] for i=1:b1) for j=1:s1)+H*w + 1e3*slack);

#for swedish market uncomment the first one
#@constraint(RobModel, noPFAD, y[3]==0)

#constraints from deterministic
@constraint(RobModel, supply_limit[i=s], x[i,j]<=PC_s[i]);
@constraint(RobModel, pretreatment[j=b], sum(a_b[j]*x[i,j] for i=1:s1)==y[j]);
@constraint(RobModel, mass_balance_blending, w==sum(y[i] for i=b));
@constraint(RobModel, demand_fullfillment, z>=D_tot);
# Fabs changed this one for less binary terms.
@constraint(RobModel, blending, (sum(Beta_b[i]*y[i] for i=1:b1))==rho*w);
@constraint(RobModel, hydro_conversion, rho*w==z);

#17

@constraint(RobModel, [ii in set_SVB],
    sum(a_b[ib] * Beta_b[ib] *
    sum(sum(B_sb[s2,ib] * D[s2,iSV] * sum((λ[s1,iSV]-μ[s1,iSV]) * Q[s1,s2] for s1 in s) for s2 in s)
    for iSV in set_SV) for ib in b)
    - η*θ[ii] >= D_tot * 100 - slack
)

#18
@constraint(RobModel, [j in s, i in set_SV],
    λ[j,i] + μ[j,i] - η * α_val[i] == 0
)


#19
@constraint(RobModel, [j in s],
    sum(sum(Q[j,j1] * (λ[j1,i] - μ[j1,i]) for j1 in s) for i in set_SV) == x[j]
)

set_optimizer_attributes(RobModel,
    "NonConvex" => 2,
    "NumericFocus" => 3,
    "OptimalityTol" => 1e-9,
    "MIPGap" => 1e-9,
    "FeasibilityTol" => 1e-9,
    "IntFeasTol" => 1e-9,
    "Presolve" => 1,
    "ScaleFlag" => 3
)

optimize!(RobModel)
@info value.(slack) == 0.0 ? "Model is feasible." : "Model is infeasible. Demand not met: $(value.(slack))"

x_rob = value.(x)
obj_rob = objective_value(RobModel)

