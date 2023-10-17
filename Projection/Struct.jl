"""
LUCA MECCA
lmecca@london.edu
July 2023
"""

#define the quadrature nodes and weights here for brevity
nodes_dict=Dict(2 => [-0.7071067811, 0.7071067811], 3 => [-1.224744871, 0, 1.224744871], 4 => [-1.650680123, -0.5246476232, 0.5246476232, 1.650680123], 5=>[-2.02018287, -0.9585724646, 0, 0.9585724646, 2.02018287],
7=>[-2.651961356, -1.673551628, -0.8162878828, 0, 0.8162878828, 1.673551628, 2.651961356], 10=>[-3.436159118, -2.5232731674, -1.756683649, -1.036610829, -0.3429013272, 0.3429013272, 1.036610829, 1.756683649, 2.5232731674, 3.436159118])
weights_dict=Dict(2 => [0.8862269254, 0.8862269254], 3 => [0.2954089751, 1.18163590, 0.2954089751], 4 => [0.08131283544, 0.8049140900, 0.8049140900, 0.08131283544], 5=>[0.01995324204, 0.3936193231, 0.9453087204, 0.3936193231, 0.01995324204],
7=>[0.0009717812450, 0.05451558281, 0.4256072526, 0.8102646175, 0.4256072526,0.05451558281, 0.0009717812450], 10=>[0.000007640432855, 0.001343645746, 0.03387439445, 0.2401386110, 0.6108626337,0.6108626337, 0.2401386110, 0.03387439445, 0.001343645746, 0.000007640432855])
    

#This file collects the constructors
#This constructor collects the steady state values
Base.@kwdef mutable struct SS
   Y::Union{Float64, Int64}
   C::Union{Float64, Int64}
   N::Union{Float64, Int64}
   WP::Union{Float64, Int64}
   X1::Union{Float64, Int64}
   X2::Union{Float64, Int64}
   Π_star::Union{Float64, Int64}
   d::Union{Float64, Int64}
end


#This constructor desribes the properties of state_variables
Base.@kwdef mutable struct state_var
    μ::Union{Float64, Int64}=0
    ρ::Union{Float64, Nothing}=nothing
    σ::Union{Float64, Int64,Nothing}=nothing
    c::Union{Float64, Int64}=0
    lb::Union{Float64, Int64, Nothing}=nothing
    ub::Union{Float64, Int64, Nothing}=nothing
    n::Int64=3
    #If no values are provided, 
    function state_var(μ::Union{Float64, Int64}=0,ρ::Union{Float64, Nothing} = nothing, σ::Union{Float64, Int64, Nothing} = nothing, c::Union{Float64, Int64} = 0, lb::Union{Float64, Int64, Nothing}=nothing, ub::Union{Float64, Int64, Nothing}=nothing, n::Int64=3)
        new(μ, ρ, σ, c, lb, ub, n)
    end
end



#This constructor contains the values that are invariant to the iterations of the algorithm
Base.@kwdef mutable struct fixed_values
    state_variables::Union{Vector{state_var}, Nothing}=nothing

    σ::Union{Int64, Float64}=1
    ϕ::Union{Int64, Float64}=5
    β::Union{Int64, Float64}=0.99
    ϵ::Union{Int64, Float64}=9
    
    α::Union{Int64, Float64}=0.25
    θ::Union{Int64, Float64}=0.75
    
    ϕ_π::Union{Int64, Float64}=1.5
    ϕ_y::Union{Int64, Float64}=0.5/4
    
    η::Union{Int64, Float64}=4

    product::String="complete"

    K::Int64=1000
    T::Int64=10000
    brns::Int64=9500
    qrts::Int64=15
    ϵ_v::Float64=0.25

    state_n::Union{Int64, Nothing}=nothing
    combination_n::Union{Int64, Nothing}=nothing
    collocation_n::Union{Int64, Nothing}=nothing

    collocation_points::Union{Matrix{Float64}, Nothing}=nothing
    Ch_bases_combinations::Union{Matrix{Float64}, Nothing}=nothing

    quadrature_n::Int64=5
    quadrature_nodes::Dict{Int64, Vector{Float64}} = nodes_dict
    quadrature_weights::Dict{Int64, Vector{Float64}} = weights_dict
   
    steady_state::Union{SS, Nothing}=nothing

    #If no values are provided, 
    function fixed_values(state_variables::Union{Vector{state_var}, Nothing}=nothing, σ::Union{Int64, Float64}=1, ϕ::Union{Int64, Float64}=5, β::Union{Int64, Float64}=0.99,
        ϵ::Union{Int64, Float64}=9, α::Union{Int64, Float64}=0.25,
        θ::Union{Int64, Float64}=0.75, ϕ_π::Union{Int64, Float64}=1.5,
        ϕ_y::Union{Int64, Float64}=0.5/4, η::Union{Int64, Float64}=4,
        product::String="complete", K::Int64=1000, T::Int64=10000, brns::Int64=9500, qrts::Int64=15, ϵ_v::Float64=0.25,
        state_n::Union{Int64, Nothing}=nothing, combination_n::Union{Int64, Nothing}=nothing,
        collocation_n::Union{Int64, Nothing}=nothing, collocation_points::Union{Matrix{Float64}, Nothing}=nothing, Ch_bases_combinations::Union{Matrix{Float64}, Nothing}=nothing,
        quadrature_n::Int64=5, quadrature_nodes::Dict{Int64, Vector{Float64}} =nodes_dict,
        quadrature_weights::Dict{Int64, Vector{Float64}} = weights_dict,
        steady_state::Union{SS, Nothing}=nothing)
        new(state_variables, σ, ϕ, β, ϵ, α, θ, ϕ_π, ϕ_y, η, product, K, T, brns, qrts, ϵ_v, state_n, combination_n,
        collocation_n, collocation_points, Ch_bases_combinations, quadrature_n, quadrature_nodes, quadrature_weights, steady_state)
    end

end


#This constructor collects the nodes and weights for the Gauss-Hermite quadrature
#Recall that weights change at each iteration
Base.@kwdef mutable struct GH_nodes_weights
    nodes_combination::Union{Vector{Matrix{Float64}}, Nothing}=nothing
    weights_combination::Union{Matrix{Float64}, Nothing}=nothing

    #If no values are provided, 
    function GH_nodes_weights(nodes_combination::Union{Vector{Matrix{Float64}}, Nothing}=nothing, weights_combination::Union{Matrix{Float64}, Nothing}=nothing)
        new(nodes_combination, weights_combination)
    end

end

#This constructor contains the unknown coefficients of the Chebyshev polynomials
Base.@kwdef mutable struct coeffs
    C_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    X1_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    Π_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing

    Π_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing

    #If no values are provided, 
    function coeffs(C_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        X1_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        Π_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        Π_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing)
        new(C_t1,X1_t1,Π_t1,Π_t)
    end
    
end

#This constructor contains the policy functions and the expectation terms
Base.@kwdef mutable struct policy_functions
    C_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    X1_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    Π_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing

    Π_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    Π_star_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing

    d_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing

    X1_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    X2_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing

    C_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    N_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    real_wage::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing

    mc_1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    mc_2::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing

    exp_X1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    exp_X2::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing
    exp_C::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing

    #If no values are provided, 
    function policy_functions(C_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        X1_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        Π_t1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        Π_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        Π_star_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        d_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        X1_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        X2_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        C_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        N_t::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        real_wage::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        mc_1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        mc_2::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        exp_X1::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        exp_X2::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing,
        exp_C::Union{Vector{Float64}, Vector{Int64}, Nothing}=nothing)
        
        new(C_t1, X1_t1, Π_t1, Π_t, Π_star_t, d_t,
        X1_t, X2_t, C_t, N_t, real_wage, mc_1, mc_2, exp_X1, exp_X2, exp_C)
    end

end