include("data_instances.jl")
include("generate_data.jl")

# Function creates a TOY problem instance with all the data. 
# Only contains data for deterministic optimization.
# 
# Data returned as struct det_ins which is then used in the optimization.
function create_toy_instance()

    # B_sb: matrix where 1 = biomass b supplied by supplier s
    B_sb=[0 0 1;0 0 1;0 0 1;0 0 1;0 1 0;0 1 0;0 1 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0];

    ins = DetInstance(
        # Sets
        3,     # b1: number of biomass  
        12,    # s1: number of suppliers 
        0:9,   # k:  decimal digits for NMDT
        1,     # p:  cardinal of l:1,2 or 3 
    
        [0.94; 0.95 ; 0.98], # a_b
        1500,                # D_tot
    
        # x_sb: amount of biomass b supplied by supplier s
        [20;30;12;10;10;20;10;15;10;10;10;10].*B_sb,            # x_sb_min
        [200;300;120;100;100;200;400;150;100;180;60;150].*B_sb, # x_sb_max
    
        # V_b: property value of biomass b
        [0.6 0.6 0.6],    # V_b_min
        [0.85 0.95 0.9],  # V_b_max
    
        # NOTE 1: rho_min 0.6 > 0.58,0.55 ? 
        # rho_sb The property value of biomass bought from supplier s
        [0.8 ;0.9; 0.78;0.64;0.68;0.85;0.7;0.58;0.7;0.65;0.55;0.6],  # rho_sb
        0.6,  # rho_min: Minimum property value of BIO-DIESEL
        0.9,  # rho_max: Maximum property value of BIO-DIESEL
    
        0.94,  # α: hydroconversion yield
    
        [90, 120, 100], # CB_b: cost of biomass b
        [50 60 40],     # PC_b: processing cost of biomass b

        # TC_sb: transportation cost of biomass b bought from supplier s 
        [100;90;80;50;100;120;100;120;120;100;80;50],  # TC_sb
        100  # HC: hydrotreatment cost 
    
    )
    return ins

end

# Function creates a GENERAL problem instance with all the data.
# Contains data for deterministic and robust optimization.
# 
# Data returned as struct 
# - det_ins which is used in deterministic optimization
# - rob_ins which is used in robust optimization
function create_deterministic_instance(suppliers, cap_factor)

    # Sets
    b1=3;         #number of biomass  
    s1=suppliers; #number of suppliers, set in main.jl
    k=0:9;        #decimal digits
    #p=accuracy;   #cardinal of l:1,2 or 3 # now set in create_linear_model(,,p)

    a_b=[0.94; 0.95 ; 0.98]; #used in overleaf 
    D_tot=1500;

    # Generate and/or update supply paramateres
    #  - PC_s: param of rob. model
    #  - B_sb: param of rob. model
    #  - D_s: param of rob. model
    #  - x_sb_min
    #  - s1
    PC_s, B_sb, D_s, x_sb_min, s1 = create_supply_and_capillarity(s1, D_tot, cap_factor)
    
    p_high=0.2 #for small suppliers
    p_low=0.05 #for big suppliers, for capillarity expansion

    # Create rho_sb with updated parameters
    # In addition, a_b, p_high, p_low params used
    #  - rho_sb The property value of biomass bought from supplier s
    rho_sb = create_qualities(p_high, p_low, PC_s, B_sb, a_b, s1); 

    # Set the remaining parameters

    x_sb_min = x_sb_min.*B_sb
    x_sb_max = PC_s.*B_sb # NOTE 2: is this correct?

    V_b_min=[0.6 0.6 0.6]; #Minimum property value of biomass b
    V_b_max=[0.85 0.95 0.9]; #Maximum property value of biomass b

    rho_min=0.6; #Minimum property value of bio-diesel
    rho_max=0.9; #Maximum property value of bio-diesel

    α=0.94; #hydroconversion yield?

    CB_b = [90, 120, 100]; #values used in overleaf
    PC_b = [50 60 40];#values used in overleaf
    T_b = 123; #transportation cost for biomass b  tons
    TC_sb=D_s*T_b; #distance*transportation costs
    HC = 100; #hydrotreatment cost 

    det_ins = DetInstance(
        b1,
        s1, # updated in capillarity expansion
        k,
        #p, set in create_linear_model(,,p)
        B_sb,

        a_b,
        D_tot,

        x_sb_min, # updated in supply and cap. exp.
        x_sb_max, # updated in supply and cap. exp.
        V_b_min,
        V_b_max,

        rho_sb, # created in create_qualities()
        rho_min,
        rho_max,
        
        α,

        CB_b,
        PC_b,
        TC_sb, # updated in supply and cap. exp.
        HC
    ) 

    # Function call to create the robust instance
    #rob_ins = create_robust_instance(s1, B_sb, data_entries, risk_v)

    #return det_ins, rob_ins;
    return det_ins
end

# Function creates the robust instance.
#
# Used in function above.
function create_robust_instance(det_ins, data_entries)

    println("calling create_robust_instance()")

    # Function call to generate robust data in generate_data.jl
    println("calling create_robust_data()")
    Km, M, D, Q = create_robust_data(det_ins, data_entries)

    # Create instance
    rob_ins = RobInstance(
        D, 
        Q,

        Km,
        M
    )

    return rob_ins
end

function create_SVC_instance(rob_ins, risk_ν)

    println("calling create_SVC_data()")
    θ, α_val, set_SV, set_SVB = create_SVC_data(rob_ins, risk_ν)

    # Create instance
    SVC_ins = SVCInstance(
        θ,
        α_val,

        set_SV,
        set_SVB
    )

    return SVC_ins
end
