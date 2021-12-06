# FILE HAS ALL DATA GENERATION FUNCTIONS

# Helper/calling function: Sets some initial data and calls functions
#  - create_supply()
#  - create_capillarity_expansion()
# defined below
#
# Returns updated supply values after both function calls
# PC_s, B_sb, D_s, x_sb_min, s1
function create_supply_and_capillarity(s1, D_tot)

    println("calling create_supply_and_capillarity()")

    # Statistics on production capacity for sampling synthetic supply data

    # PFAD
    sdPFAD=11.5; #palm oil analytics database
    avgPFAD=76;
    sdMINPFAD=9; #calculated from toy numbers in overleaf
    avgMINPFAD=18; #The minimun values are taken from overleaf, not properly sourced 

    # Animal 
    # European comission trade database, EU comission report of emissions of animal fats
    # used for biodiesel fats
    sdAF=100; 
    avgAF=1500;
    sdMINAF=13;
    avgMINAF=5;

    # Cooking oils
    sdUCO=100; #European comission trade database and UK biofuel statistics
    avgUCO=500;
    sdMINWCO=11;
    avgMINWCO=2;

    # Generate supply parameters
    #  - PC_s
    #  - B_sb
    #  - D_s
    #  - x_sb_min
    # with function create_supply()
    PC_s, B_sb, D_s, x_sb_min = create_supply(s1, D_tot,
                                              sdPFAD, sdAF, sdUCO,
                                              avgPFAD, avgAF, avgUCO,
                                              sdMINPFAD, avgMINPFAD,
                                              sdMINAF, avgMINAF,
                                              sdMINWCO, avgMINWCO)

    capillarity_factor=2; #give as input, test different levels

    # Implement capillarity expansion for CO suppliers
    # Updates supply parameters
    #  - PC_s
    #  - B_sb
    #  - D_s
    #  - x_sb_min
    #  - s1
    # with function create_capillarity_expansion()
    PC_s, B_sb, D_s, x_sb_min, s1 = create_capillarity_expansion(capillarity_factor, PC_s, B_sb, D_s, x_sb_min, s1)

    return PC_s, B_sb, D_s, x_sb_min, s1
end

# Function creates synthetic supply data drawn from distributions
# with means and standard deviations from above
#
# s = number of suppliers
# D_tot = total demand
#
# Function returns:
#  - PC_s: production capacities of suppliers
#  - B_sb: matrix indicating which suppliers produce which biomass
#  - D_s:  distances to suppliers
#  - x_sb_min: minimum amount of biomass b expected from supplier s
function create_supply(s, D_tot, sdPFAD, sdAF, sdUCO, avgPFAD, avgAF, avgUCO, sdMINPFAD, avgMINPFAD, sdMINAF, avgMINAF,sdMINWCO, avgMINWCO)

    # Set nr of suppliers per biomass
    PFAD=rand(DiscreteUniform(1, s÷3)) 
    WCO=rand(DiscreteUniform(s÷3, s÷2))
    AF=s-PFAD-WCO

    # Set distances, production capacities and _ to suppliers
    Ds_test=zeros(Float64, s)
    PC_test=zeros(Float64, s)
    x_sb_min=zeros(Float64, s) 

    ######################### creation of B_sb
    
    # B_sb_test is a sx3 matrix where [i,j] is 1 if 
    # supplier i produces biomass j
    B_sb_test=zeros(Int64, s, 3) # sx3 matrix of zeroes
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
    
    # Distances to biomass suppliers are sampled from 
    # different normal distributions
    for i in 1:PFAD
        dist_sample=rand(Normal(1100, 150)) #*10 THESE ARE ONE WAY, move numbers out of the function?
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
    while( sum(PC_test) < D_tot ) # ensure total PC meets total demand
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


    ######################### creation of x_sb_min
    for i in 1:PFAD
        min_sample=rand(Normal(avgMINPFAD, sdMINPFAD)) 
        if (min_sample<0)
            min_sample=0
        end
        x_sb_min[i]=min_sample
    end
    
    for k in PFAD+1:PFAD+AF
        min_sample=rand(Normal(avgMINAF, sdMINAF)) 
        if (min_sample<0)
            min_sample=0
        end
        x_sb_min[k]=min_sample
    end
     
    for j in PFAD+AF+1:s
        sample=rand(Normal(avgMINWCO, sdMINWCO)) 
        if (sample<0)
            sample=0
        end
        x_sb_min[j]=sample
    end
     
    
    return round.(PC_test, digits=2), round.(B_sb_test, digits=2), round.(Ds_test, digits=2), round.(x_sb_min, digits=2) 
end 

# Function increases the nr of CO suppliers by 'capillarity_factor'.
# The sum of PC, x_sb_min remain. 
#
# Parameters from function above + cap.fact., additionally:
#  - r_low, r_high: radiuses used for updated distances
#
# Returns:
#  - newB_sb: new CO suppliers added
#  - newPC_s: PC_s divided to new CO suppliers
#  - newD_s: new distances, added rand radius to each
#  - new_s1: updated nr of suppliers
#  - new_s: 1:new_s1
#  - new_x_sb_mins: x_sb_mins divided to new CO suppliers
#
# CHANGE: splitting randomly not to even parts

function create_capillarity_expansion(capillarity_factor, PC_s, B_sb, D_s, x_sb_min, s1, r_low=50, r_high=100)

    ####################### Update nr of CO suppliers from B_sb into newB_sb
    UCO_suppliers=sum(B_sb[:,1]) # nr of UCO suppliers
    UCO_suppliers=Int(UCO_suppliers)
    
    #how many more UCO will be, all suppliers divided to #capillarity_factor many suppliers
    UCO_new=capillarity_factor*UCO_suppliers # updated nr of CO suppliers
    newB_sb=B_sb[1:(s1-UCO_suppliers),:] # B_sb without CO suppliers

    matrix=zeros(Int8, UCO_new , 3) # matrix for new CO suppliers
    matrix[:,1].=1
    newB_sb=vcat(newB_sb,matrix) # add updated nr of CO suppliers back to newB_sb
    #######################
    
    #now go through PC_s of UCO suppliers and explode them to little suppliers
    # Divide original PC to larger number of CO suppliers, total PC remains the same
    copy=PC_s[s1-UCO_suppliers+1:s1] # PC_s of CO suppliers
    newPC_s=PC_s[1:s1-UCO_suppliers] # PC_s other suppliers
    for j in 1:length(copy) # for each CO supplier PC
        for k in 1:capillarity_factor
            value=copy[j]/capillarity_factor # randomly split the capacity, five points
            append!(newPC_s, value) 
        end
    end

    #now go through x_sb_mins of UCO suppliers and explode them to little suppliers
    # Divide original x_sb_mins to larger number of CO suppliers, total x_sb_mins remains the same
    copyMIN=x_sb_min[s1-UCO_suppliers+1:s1] # x_sb_min of CO suppliers
    newMINS=x_sb_min[1:s1-UCO_suppliers] # x_sb_min of other suppliers
    for j in 1:length(copyMIN)
        for k in 1:capillarity_factor
            valueMIN=copyMIN[j]/capillarity_factor # randomdy split the capacity, five points
            append!(newMINS, valueMIN)
        end
    end

    # Add some radius to original distances
    DsCopy=D_s[s1-UCO_suppliers+1:s1] #only uco distances
    newD_s=D_s[1:s1-UCO_suppliers]    #no uco distances
    
    for a in 1:length(DsCopy) # for old distances
        for b in 1:capillarity_factor
            radius=rand(DiscreteUniform(r_low, r_high))
            val=DsCopy[a]+radius # update distance with old dist + random radius
            append!(newD_s, val)
        end
    end

    new_s1=s1-UCO_suppliers+UCO_new
    return newPC_s, newB_sb, newD_s, newMINS, new_s1 
end



# Function counts the property value (rho) of biomass bought from supplier s
#  - property value: ratio describing how much biodiesel we get from initial biomass 
#  - the property value is sampled from a norm distribution 
#
#  - mean: the preprocessing yield of biomass
#  - std: the preprocessing yield of biomass * prop. value of BIODIESEL
#      - prop. value of BIODIESEL: describes the yield of blending
#      - the values are between p_low and p_high
#      - and assigned based on size (PC) of supplier
#
# Parameters: PC_s, B_sb, s1, additionally
#  - p_high: best property value
#  - p_low: worst property value
#  - a_b: preprocessing yield of biomass types [0.94; 0.95; 0.98]
#
# Returns:
#  - q_s: rho_s: prop value of blend. treated biomass
#
function create_qualities(p_high, p_low, PC_s,B_sb,a_b,s1)

    println("calling create_qualities()")

    s = 1:s1

    interval=(p_high-p_low)/(s1-1) 
    q_s=zeros(Float64,s1) #quality dimension

    for i in s # def on last func call, for each i in nr of supplier

        # Find sorted index of this PC
        j=1
        while PC_s[i] != (sort(PC_s))[j]
            j=j+1
        end
        
        index = j

        # B_sb*a_b = preprocessing yield of biomass b from supplier s
        q_s[i] = rand(Normal((B_sb*a_b)[i], (p_high - (index-1)*interval)*(B_sb*a_b)[i] )) 
        
        if q_s[i]>1
            q_s[i]=1
        end
    end
    return q_s
end

# Function generates all the data for the robust data instance.
# Takes as input
#  - s1: number of suppliers
#  - B_sb: binary matrix mapping suppliers to biomasses
#  - data_entries: number of entries to use in SVC model
#
# Returns the robust data
function create_robust_data(s1, B_sb, data_entries)

    s = 1:s1

    # Data creation for robust model

    #gamma_s: uncertainty factor for each supplier [0,1]
    #gamma matrix Sx100 historical data for all suppliers
    ent = data_entries;
    D = zeros(Float64, length(s), ent); #saving gamma data here

    # Sample ent nr of gammas for all suppliers
    # Distr. varies on type of biomass (based on historical data)
    for i in s
        # UCO 
        if B_sb[i,1]==1
            Gamma_b=round.(rand(Normal(75, 5), ent), digits=2) #WCO
            D[i,:]=Gamma_b

        # AF
        elseif B_sb[i,2]==1
            Gamma_b=round.(rand(Normal(95, 5), ent), digits=2)   #AF
            D[i,:]=Gamma_b

        # PFAD
        else B_sb[i,3]==1
            Gamma_b=round.(rand(Normal(85, 5), ent), digits=2)   #PFAD
            D[i,:]=Gamma_b
        end
    end

    # Limit values to [0.00, 100.00]
    for i in s
        for j in 1:ent
            D[i,j] = (D[i,j]>100.00) ? 100.00 : D[i,j]
            D[i,j] = (D[i,j]<0.00) ? 0.00 : D[i,j]
        end
    end
    # D is now gamma_s s x entities

    Σ = cov(D, corrected=true, dims=2);
    Q = round.(Σ^(-0.5), digits = 6);

    println("Type of Q before: $(typeof(Q)))")

    println()
    println("Diagonal elems of Σ[1:10] before")
    Σ_diag = Σ[diagind(Σ)]
    for elem in Σ_diag[1:10] 
        println(elem)
    end

    # Q matrix sometimes gets complex numbers
    # Make small modifications to the diagonal of Σ
    # to avoid this
    while typeof(Q) != Matrix{Float64}
        Σ[diagind(Σ)] .+= 1e-5
        Q = round.(Σ^(-0.5), digits = 6);
    end

    println()
    println("Type of Q after: $(typeof(Q)))")

    println("Diagonal elems of Σ[1:10] after")
    Σ_diag = Σ[diagind(Σ)]
    for elem in Σ_diag[1:10] 
        println(elem)
    end

    # Step 2: form the Kernel matrix

    # Min and max gamma for each supplier
    uU = maximum(D, dims=2); #dims=2 since data is per column
    uL = minimum(D, dims=2);

    N = s1; #N is the same as the number of suppliers
    J = 1:N ;

    # Generic function K
    #  - sums the diff. of gamma (sample) entries u,v to gamma min max for each supplier
    #  - u,v: parameter vectors (gamma entry vectors for suppliers)
    #  - uU[j] - uL[j]: difference of max and min gamma for supplier j (0.95-0.60)
    #  - abs(u[v] - v[j]): diff. of gamma entry u and v for supplier j abs(0.83-0.86)
    K(u,v) = sum(uU[j] - uL[j] - abs(u[j] - v[j]) for j in J);


    # Construct ent x ent kernel matrix Km with K(u,v)

    # The kernel trick K sweeps the data to be calculated.
    M = ent; #data entities
    Km = zeros(M, M);
    for i = 1:M
        for i1 = 1:M
            # D[:,i], D[:,i1] = sx1 vector of gammas of each supplier 
            Km[i,i1] = K(D[:,i], D[:,i1])
        end
    end


    # Step 3: derive the support vector clustering (SVC) model

    ν = 0.05;   # regularisation parameter

    SVC_Model = Model(Gurobi.Optimizer);
    #set_optimizer_attributes(SVC_Model, "mehrotra_algorithm" => "yes");

    # Weight variable alpha
    @variable(SVC_Model, α[i in 1:M] >= 0);

    # 37) NOTE 3: a_j in report?
    @objective(SVC_Model, Min,
        sum(α[i] * α[j] * Km[i,j] for i in 1:M for j in 1:M)
        - sum(α[i] * Km[i,i] for i in 1:M)
    );

    @constraint(SVC_Model, sum(α[i] for i in 1:M) == 1); # 38)
    @constraint(SVC_Model, [i in 1:M], α[i] * (M * ν) <= 1); # 39)

    println(" ----- SVC model ----- ")
    optimize!(SVC_Model)
    println(" ----- End of of SVC model ----- ")

    # Format optimized a
    α_val = round.(value.(α), digits = 6);
    for i in 1:length(α_val)
        α_val[i] < 1e-6 ? α_val[i] = 0.0 : α_val[i]
    end

    println()
    println("Values of a_val[1:10]:")
    for elem in α_val[1:10]
        println(elem)
    end

    # Filtering the support vectors and the boundary support vectors
    #  - SV: values at lower boundary 0, the rest 1 
    #  - SVB: values at upper and lower boundary 0, the rest 1
    ϵ = 1e-12 # tolerance for zero: helps finding boundary points correctly.
    SV = [α_val[i] > 0 + ϵ ? 1 : 0 for i in 1:M]
    SVB = [α_val[i] > 0 + ϵ  && α_val[i] < 1/(M * ν) - ϵ ? 1 : 0 for i in 1:M]

    # Forming the sets with (boundary) support vectors
    # Collect indexes of 1s
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

    θ = Dict(); # NOTE 5: what is this

    for i1 in set_SVB # = for indexes of 1 in SVB
        θ[i1] = sum(α_val[i] * norm(Q * (D[:,i1] - D[:,i]), 1) for i in set_SV)
    end

    return D, Q, θ, α_val, set_SV, set_SVB

end
