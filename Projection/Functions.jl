"""
LUCA MECCA
lmecca@london.edu
July 2023
"""

#This file contains all the time-intensive functions specific to this project 
#This functions solves for the steady state of output
function compute_SS(constant_values::fixed_values):SS
    ϕ=constant_values.ϕ
    ϵ=constant_values.ϵ
    θ=constant_values.θ
    α=constant_values.α
    β=constant_values.β
    σ=constant_values.σ

    function g!(G, ζ)
        Y,C, N, WP, X1, X2, Π_star,d=ζ[1],ζ[2],ζ[3],ζ[4],ζ[5],ζ[6],ζ[7],ζ[8]

        G[1]=C^σ*N^ϕ-WP
        G[2]=ϵ*X1-(ϵ-1)*X2*Π_star
        G[3]=X1-WP/(N^(-α)*(1-α))*Y/C^σ-θ*β*X1
        G[4]=X2-Y/C^σ-β*θ*X2
        G[5]=1-θ-(1-θ)*Π_star^(1-ϵ)
        G[6]=d-θ*d-(1-θ)*Π_star^(-ϵ/(1-α))
        G[7]=Y-C
        G[8]=N-Y^(1/(1-α))*d

    end

    ζ=ones(8)
    sol=nlsolve(g!, ζ, method = :newton)
    steady_state=SS(sol.zero[1], sol.zero[2], sol.zero[3], sol.zero[4], sol.zero[5], sol.zero[6], sol.zero[7], sol.zero[8])
    return steady_state
end



#Inflation for the fraction of firms that are allowed to change prices
function Π_star(constant_values::fixed_values, policy_Π::Vector{Float64})::Vector{Float64}

    θ=constant_values.θ
    ϵ=constant_values.ϵ

    P_change_star=((1 .-θ.*policy_Π.^(ϵ-1))./(1-θ)).^(1/(1-ϵ))
    return P_change_star
end

#Second auxiliary variable
function x2_update(constant_values::fixed_values, policy_X1::Vector{Float64}, policy_Π_star::Vector{Float64})::Vector{Float64}
    ϵ=constant_values.ϵ
    α=constant_values.α
    x2=ϵ/(ϵ-1).*policy_X1./policy_Π_star.^(1+ϵ*α/(1-α))
    return x2
end

#Second auxiliary variable
function x1_update(constant_values::fixed_values, policy_X2::Vector{Float64}, policy_Π_star::Vector{Float64})::Vector{Float64}
    ϵ=constant_values.ϵ
    α=constant_values.α
    x1=(ϵ-1)/ϵ.*policy_X2.*policy_Π_star.^(1+ϵ*α/(1-α))
    return x1
end

#Law of motion of price dispersion d_t
function law_motion_d(constant_values::fixed_values, policy_f::policy_functions)::Vector{Float64}
    θ=constant_values.θ
    ϵ=constant_values.ϵ
    α=constant_values.α
    collocation_points=constant_values.collocation_points
    policy_Π=policy_f.Π_t
    policy_Π_star=policy_f.Π_star_t
    state=constant_values.state_variables[4]

    collocation_d_X=Ch_to_state(state, (collocation_points)[4,:])

    d_update=θ.*policy_Π.^(ϵ/(1-α)).*collocation_d_X .+ (1-θ).*policy_Π_star.^(-ϵ/(1-α))
    return d_update    
end

#This function computes the second auxiliary variable x_{2,t} from its intertemporal equation
function x2_int(constant_values::fixed_values, policy_f::policy_functions)::Vector{Float64}
    β=constant_values.β
    θ=constant_values.θ
    σ=constant_values.σ
    state_Z=constant_values.state_variables[2]
    Z_t=exp.(Ch_to_state(state_Z,constant_values.collocation_points[2,:])) #Preference shock Z_t
    
    C_t=policy_f.C_t
    exp_term=policy_f.exp_X2

    policy_x2= C_t.^(1-σ).*Z_t .+ β*θ .* exp_term
    return policy_x2
end

#This function computes consumption from the Euler Equation
function C_EE(constant_values::fixed_values, policy_f::policy_functions)::Vector{Float64}
    ϕ_π=constant_values.ϕ_π
    ϕ_y=constant_values.ϕ_y
    σ=constant_values.σ
    state_Z=constant_values.state_variables[2]
    state_V=constant_values.state_variables[3]
    collocation_points=constant_values.collocation_points
    policy_Π=policy_f.Π_t
    exp_term=policy_f.exp_C
  
    constant_values_SS=fixed_values()
    constant_values_SS.state_variables, constant_values_SS.state_n=all_state_variables, lastindex(all_state_variables)
    constant_values_SS.collocation_points=[0,0,0,state_to_Ch(constant_values_SS.state_variables[4],1)][:,:] #steady state
    constant_values_SS.collocation_n, constant_values_SS.combination_n=size(constant_values_SS.collocation_points,2), constant_values.combination_n
    constant_values_SS.Ch_bases_combinations=Ch_bases_function(constant_values_SS) #bases functions

    Y_ss=exp.(Ch_pol(constant_values_SS, coefficients.C_t1))[1]

    Z_t=exp.(Ch_to_state(state_Z,collocation_points[2,:])) #Preference shock Z_t
    V_t=exp.(Ch_to_state(state_V,collocation_points[3,:])) #Monetary policy shock V_t

    C_t=(exp_term./Y_ss^ϕ_y.*policy_Π.^ϕ_π.*V_t./Z_t).^(-1/(σ+ϕ_y))

    return C_t
end

#This functions computs the level of employment from the level fo aggregate supply
function aggregate_supply(constant_values::fixed_values, policy_f::policy_functions)::Vector{Float64}
    α=constant_values.α
    state=constant_values.state_variables[1]
    policy_Y=policy_f.C_t
    policy_d=policy_f.d_t
    collocation_points=constant_values.collocation_points

    A_t=exp.(Ch_to_state(state,collocation_points[1,:])) #Technology shock A_t
    
    N_t=(policy_Y./A_t).^(1/(1-α)).*policy_d
    return N_t
end

#This function computes the real wage from the labour-leisure choice
function LL_choice(constant_values::fixed_values, policy_f::policy_functions)::Vector{Float64}
    ϕ=constant_values.ϕ
    policy_C=policy_f.C_t
    policy_N=policy_f.N_t
    σ=constant_values.σ

    real_wage=policy_C.^σ.*policy_N.^ϕ
    return real_wage
end

#This funcitons computes the real marginal cost from the policy functions for real wage and employment
function mc_hat(constant_values::fixed_values, policy_f::policy_functions)::Vector{Float64}
    α=constant_values.α
    state=constant_values.state_variables[1]
    policy_N=policy_f.N_t
    policy_d=policy_f.d_t
    policy_real_wage=policy_f.real_wage
    collocation_points=constant_values.collocation_points

    A_t=exp.(Ch_to_state(state,collocation_points[1,:])) #Technology shock A_t
    
    marginal_cost=policy_real_wage./(A_t.*policy_N.^(-α).*(1-α).*policy_d.^(α-1))

    return marginal_cost
end

#This functions computes the real marginal cost from the equilibrium condition of the first auxiliary variable
function mc(constant_values::fixed_values, policy_f::policy_functions)::Vector{Float64}
    θ=constant_values.θ
    β=constant_values.β
    σ=constant_values.σ
    C_t=policy_f.C_t
    policy_X1=policy_f.X1_t
    policy_d=policy_f.d_t
    exp_term=policy_f.exp_X1
    state_Z=constant_values.state_variables[2]
    Z_t=exp.(Ch_to_state(state_Z,constant_values.collocation_points[2,:])) #Preference shock Z_t

    marginal_cost=(policy_X1.-β*θ.*exp_term)./(Z_t.*C_t.^(1-σ).*policy_d)
    return marginal_cost
end


#Function Tau discretizes AR(1) processes into discrete grids using the approach proposed in Tauchen (1986)
#It also computes the corresponding transition matrix
#Inputs:
#ρ: autoregressive coefficient.
#σ: volatility of the shock.
#K: number of grid points. Default value is 9.
#λ: MC approx parameter. Default value is 3 (as suggested in Tauchen).
#interc: fix part. Default value is 0.

function Tau(ρ::Float64, σ::Union{Float64, Int64}, K::Int64=5, λ::Number=3, interc::Number=0)::Tuple{Vector{Float64}, Matrix{Float64}}
    #Initial point
    z_1=-λ*σ/sqrt(1-ρ^2)
    #size of the step
    step=-2*z_1/(K-1)
    #create the grid
    grid=Array{Number}(undef,K,1) 

    grid[1]=z_1  
    #complete the grid
    for i in 2:K
        grid[i]=z_1+(i-1)*step
    end 

    #add the constant part
    grid=grid .+ interc

    #Now create the transition matrix
    #T_ij indicates the probability of moving from state j to state i
    #For example the first column of the matrix includes the probability of moving from
    #the first state to itself and all the other states
    #The first row includes the probability of going from all the states to the first state.
    T_matrix = Matrix{Float64}(undef,K,K) 
    #loop over each column of the transition matrix
    for col in 1:K
        #probability of moving from any of the states to state 1 (i.e first row of the transition matrix)
        T_matrix[1,col]=cdf(Normal(),(grid[1]-interc-grid[col]*ρ)/σ+step/2/σ)
        #Compute the transition probability of in-betweeen states 
        for row in 2:(K-1)
            T_matrix[row,col]=cdf(Normal(),(grid[row]-interc-grid[col]*ρ)/σ+step/2/σ)-cdf(Normal(),(grid[row]-interc-grid[col]*ρ)/σ-step/2/σ)
        #probability of moving from any of the states to state m (i.e last row of the transition matrix)
        T_matrix[K,col]=max(1-sum(T_matrix[1:K-1,col]),0)
        end
    end
    return vec(grid), T_matrix
end


#This functions builds a set of unique points at which the policy functions (chebyshev polynomials) can be evaluated
#The points are returned in the Chebyshev donimnion [-1,1]
function evaluation_points(constant_values::fixed_values)::Matrix{Float64}
    state_variables=constant_values.state_variables
    state_n=constant_values.state_n

    vector_grid=Vector{Vector{Float64}}(undef, state_n) #container for the grids

    for i in 1:state_n
        state=state_variables[i]
        #We set the endogenous state variable price dispersion to its steady state level of 1 because it is assumed to be so in the log-linear solution
        if state==d_lag
            vector_grid[i]=state_to_Ch(state,ones(length(vector_grid[1])))
        else
            vector_grid[i]=state_to_Ch(state,Tau(constant_values.state_variables[i].ρ, constant_values.state_variables[i].σ)[1])
        end
    end

    #compute all the possible combinations and remove duplicates
    tensor_combinations=[]
    tensor_zeros(vector_grid, [], tensor_combinations)
    matrix = hcat(tensor_combinations...)
    unique_cols = Set{Vector}(eachcol(matrix))  # remove duplicates
    unique_matrix = hcat(unique_cols...)

    return unique_matrix

end
