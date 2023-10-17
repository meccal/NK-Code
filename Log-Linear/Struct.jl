"""
LUCA MECCA
lmecca@london.edu
July 2023
"""
#This file collects the data structures

#This data structure contains the calibrated parameters of the model
struct fixed_parameters
    σ::Union{Int64, Float64}
    ϕ::Union{Int64, Float64}
    β::Union{Int64, Float64}
    ϵ::Union{Int64, Float64}
    
    ρ_z::Float64
    σ_z::Union{Int64, Float64}
    n_z::Int64

    α::Union{Int64, Float64}
    θ::Union{Int64, Float64}
    
    ρ_a::Float64
    σ_a::Union{Int64, Float64}
    n_a::Int64

    ϕ_π::Union{Int64, Float64}
    ϕ_y::Union{Int64, Float64}
    
    ρ_v::Float64
    σ_v::Union{Int64, Float64}
    n_v::Int64

    η::Union{Int64, Float64}

    qrts::Int64
    ϵ_v::Union{Int64, Float64}
    ϵ_z::Union{Int64, Float64}
    ϵ_a::Union{Int64, Float64}
end

#This data structure contains the current state of the economy, as summarized by the three exogenous state variables
struct state_vars
   a_t::Union{Vector{Int64}, Vector{Float64}}
   z_t::Union{Vector{Int64}, Vector{Float64}}
   v_t::Union{Vector{Int64}, Vector{Float64}}  
end

#This data structure contains the log-linear coefficients
struct log_lin_coeff
    A1Y::Union{Int64, Float64}
    A2Y::Union{Int64, Float64}
    A3Y::Union{Int64, Float64}
    A1π::Union{Int64, Float64}
    A2π::Union{Int64, Float64}
    A3π::Union{Int64, Float64}
 end