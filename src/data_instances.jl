# Structure Instance contains all data of a problem
# to be optimized
mutable struct DetInstance 
    b1::Int64
    s1::Int64
    k::UnitRange{Int64}
    p::Int64

    a_b::Vector{Float64}
    D_tot::Int64

    x_sb_min::Matrix{Float64} # Int64 in toy
    x_sb_max::Matrix{Float64} # Int64 in toy

    V_b_min::Matrix{Float64}
    V_b_max::Matrix{Float64}

    rho_sb::Vector{Float64}
    rho_min::Float64
    rho_max::Float64
    
    α::Float64

    CB_b::Vector{Float64} # Int64 in toy
    PC_b::Matrix{Float64} # Int64 in toy
    TC_sb::Vector{Float64} # Int64 in toy
    HC::Int64
end

# NOT READY
mutable struct RobInstance
    B_sb::Matrix{Float64}
    PC_s::Vector{Float64}
    D_s::Vector{Float64}

    D::Matrix{Float64}
    Q::Matrix{Float64}
    θ::Dict{Any, Any}
    a_val::Vector{Float64}

    set_SV::Vector{Any}
    set_SVB::Vector{Any}

    #=
    C_b
    T_b
    P_b

    H
    =#
end
