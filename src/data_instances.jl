# Structure Instance contains all data of a problem
# to be optimized
mutable struct DetInstance 
    b1::Int64
    s1::Int64
    k::UnitRange{Int64}
    p::Int64

    a_b::Vector{Float64}
    D_tot::Int64

    x_sb_min::Matrix{Int64}
    x_sb_max::Matrix{Int64}

    V_b_min::Matrix{Float64}
    V_b_max::Matrix{Float64}

    rho_sb::Vector{Float64}
    rho_min::Float64
    rho_max::Float64
    
    α::Float64

    CB_b::Vector{Int64}
    PC_b::Matrix{Int64}
    TC_sb::Vector{Int64}
    HC::Int64
end

# NOT READY
#=
mutable struct RobInstance
    b1::Int64
    s1::Int64

    a_b::Vector{Float64}
    D_tot::Int64

    B_sb
    PC_s
    D_s
    D_tot

    C_b
    P_b

    set_SV
    set_SVB

    D
    Q
    θ
    a_val

    H
end
=#