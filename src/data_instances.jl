# Structure DetInstance contains all the
# deterministic data of a problem to be optimized
mutable struct ToyInstance 
    b1::Int64
    s1::Int64
    k::UnitRange{Int64}

    B_sb::Matrix{Float64}
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

mutable struct DetInstance 
    b1::Int64
    s1::Int64
    k::UnitRange{Int64}

    B_sb::Matrix{Float64}
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

# Structure RobInstance contains all the
# robust data of a problem to be optimized
mutable struct RobInstance
    D::Matrix{Float64} # uncertainty matrix
    Q::Matrix{Float64} # covariance matrix

    Km
    M
end

mutable struct SVCInstance
    θ::Dict{Any, Any}
    α_val::Vector{Float64}

    set_SV::Vector{Any}
    set_SVB::Vector{Any}
end
