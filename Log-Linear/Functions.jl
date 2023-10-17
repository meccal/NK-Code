"""
LUCA MECCA
lmecca@london.edu
July 2023
"""

#OPTIMAL MARKUP
#This function computes the optimal markup
#ϵ is the elasticity of substitution between consumption goods
function μ(parameters::fixed_parameters)::Float64
    ϵ=parameters.ϵ

    opt_markup=log(ϵ/(ϵ-1))
    return opt_markup
end


#RESPONSE OF INFLATION TO THE OUTPUT GAP
function κ(parameters::fixed_parameters)::Float64
    θ=parameters.θ
    α=parameters.α
    ϵ=parameters.ϵ
    β=parameters.β
    σ=parameters.σ
    ϕ=parameters.ϕ

    kappa=(1-θ)*(1-β*θ)/θ*((1-α)/(1-α+α*ϵ))*(σ+(ϕ+α)/(1-α))

    return kappa

end


function Ω(parameters::fixed_parameters)::Float64
    σ=parameters.σ
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π  

    Omega=1/(σ+ϕ_y+κ(parameters)*ϕ_π)
    return Omega
end


function ψ_ya(parameters::fixed_parameters)::Float64
    σ=parameters.σ
    ϕ=parameters.ϕ
    α=parameters.α  

    psi_ya=(1+ϕ)/(σ*(1-α)+ϕ+α)
    return psi_ya
end

#NATURAL RATE OF INTEREST
#The natural rate of intrest depends on both the technology and preference shocks
function r_natural(parameters::fixed_parameters, state::state_vars)::Vector{Float64}
    β=parameters.β
    σ=parameters.σ
    ρ_a=parameters.ρ_a
    ϕ=parameters.ϕ
    α=parameters.α
    ρ_z=parameters.ρ_z

    a_t=state.a_t
    z_t=state.z_t

    r_t_n=-log(β)-σ*(1-ρ_a)*(1+ϕ)/(σ*(1-α)+ϕ+α).*a_t .+(1-ρ_z).*z_t
    return r_t_n
end


#NATURAL LEVEL OF OUTPUT
#The natural level of output depends only on the technology shock, not on monetary policy and preference shocks
function y_natural(parameters::fixed_parameters, state::state_vars)::Vector{Float64}
    α=parameters.α
    σ=parameters.σ
    ϕ=parameters.ϕ

    a_t=state.a_t

    y_t_n=(1+ϕ)/(σ*(1-α)+ϕ+α).*a_t .- (1-α)*(μ(parameters)-log(1-α))/(σ*(1-α)+ϕ+α)
    return y_t_n
end


#EMPLOYMENT
#y_t is output
#a_t is the technology shock
function empl(parameters::fixed_parameters, state::state_vars, y_t::Vector{Float64})::Vector{Float64}
    α=parameters.α
    a_t=state.a_t

    n_t=(y_t-a_t)./(1-α)
    return n_t
end


#REAL WAGE
#y_t is output
#a_t is the technology shock
function real_wage(parameters::fixed_parameters, y_t::Vector{Float64}, n_t::Vector{Float64})::Vector{Float64}
    σ=parameters.σ
    ϕ=parameters.ϕ

    w_t_real=σ.*y_t+ϕ.*n_t
    return w_t_real
end

##############################################################################################
################################# IMPULSE RESPONSE FUNCTIONS #################################
##############################################################################################

########################################## MP SHOCK ##########################################

#Effect of monetary policy shock on otuput (gap)
#annual:input "Yes" if you want the IRF to be annualized, "No" otherwise
function IRF_MP_Y(parameters::fixed_parameters, state::state_vars, annual::String="No")::Vector{Float64}
    β=parameters.β
    ρ_v=parameters.ρ_v
    σ=parameters.σ
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π
   
    MP_shocks=state.v_t

    if annual=="No"
        Y_MP=(-(1-β*ρ_v)/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(parameters)*(ϕ_π-ρ_v))).*MP_shocks
    elseif annual=="Yes"
        Y_MP=(-(1-β*ρ_v)/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(parameters)*(ϕ_π-ρ_v))).*MP_shocks.*4
    else
        return(error("The parameter annual should be either Yes or No"))
    end

    return Y_MP
end


#Effect of monetary policy shock on inflation
#annual:input "Yes" if you want the IRF to be annualized, "No" otherwise
function IRF_MP_π(parameters::fixed_parameters, state::state_vars, annual::String="No")::Vector{Float64}
    β=parameters.β
    ρ_v=parameters.ρ_v
    σ=parameters.σ
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π

    MP_shocks=state.v_t
   
    if annual=="No"
        π_MP=(-κ(parameters)/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(parameters)*(ϕ_π-ρ_v))).*MP_shocks
    elseif annual=="Yes"
        π_MP=(-κ(parameters)/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(parameters)*(ϕ_π-ρ_v))).*MP_shocks.*4
    else
        return(error("The parameter annual should be either Yes or No"))
    end

    return π_MP
end


#Effect of monetary policy shock on the real interest rate
#annual:input "Yes" if you want the IRF to be annualized, "No" otherwise
function IRF_MP_r(parameters::fixed_parameters, state::state_vars, annual::String="No")::Vector{Float64}
    β=parameters.β
    ρ_v=parameters.ρ_v
    σ=parameters.σ
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π

    MP_shocks=state.v_t
   
    if annual=="No"
        r_MP=(σ*(1-ρ_v)*(1-β*ρ_v)/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(parameters)*(ϕ_π-ρ_v))).*MP_shocks
    elseif annual=="Yes"
        r_MP=(σ*(1-ρ_v)*(1-β*ρ_v)/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(parameters)*(ϕ_π-ρ_v))).*MP_shocks.*4
    else
        return(error("The parameter annual should be either Yes or No"))
    end

    return r_MP
end


#Effect of monetary policy shock on the nominal interest rate
#annual: input "Yes" if you want the IRF to be annualized, "No" otherwise
function IRF_MP_i(parameters::fixed_parameters, state::state_vars, annual::String="No")::Vector{Float64}
    β=parameters.β
    ρ_v=parameters.ρ_v
    σ=parameters.σ
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π

    MP_shocks=state.v_t
   
    if annual=="No"
        i_MP=((σ*(1-ρ_v)*(1-β*ρ_v)-ρ_v*κ(parameters))/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(parameters)*(ϕ_π-ρ_v))).*MP_shocks
    elseif annual=="Yes"
        i_MP=((σ*(1-ρ_v)*(1-β*ρ_v)-ρ_v*κ(parameters))/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(parameters)*(ϕ_π-ρ_v))).*MP_shocks.*4
    else
        return(error("The parameter annual should be either Yes or No"))
    end

    return i_MP
end


#Effect of monetary policy shock on money supply
#annual: input "Yes" if you want the IRF to be annualized, "No" otherwise
#p_t is the price level
function IRF_MP_Ms(parameters::fixed_parameters, state::state_vars, p_t::Vector{Float64}, annual::String="No")::Vector{Float64}
    β=parameters.β
    ρ_v=parameters.ρ_v
    σ=parameters.σ
    ϕ_y=parameters.ϕ_y
    η=parameters.η
    ϕ_π=parameters.ϕ_π

    MP_shocks=state.v_t
    lag_price=lag(p_t,1)
    lag_price[1]=0
   
    if annual=="No"
        Ms_MP=lag_price.-(((1-β*ρ_v)*(1+η*σ*(1-ρ_v))+(1-η*ρ_v)*κ(parameters))/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(parameters)*(ϕ_π-ρ_v))).*MP_shocks
    elseif annual=="Yes"
        Ms_MP=lag_price-(((1-β*ρ_v)*(1+η*σ*(1-ρ_v))+(1-η*ρ_v)*κ(parameters))/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(parameters)*(ϕ_π-ρ_v))).*MP_shocks.*4
    else
        return(error("The parameter annual should be either Yes or No"))
    end

    return Ms_MP
end


##############################################################################################
################################## UNDETERMINED COEFFICIENTS #################################
##############################################################################################
#The following functions compute the coefficients of the following guesses for the policy functions:
#̃y_t=A1_y a_t + A2_y z_t + A3_y v_t
#π_t=A1_π a_t + A2_π z_t + A3_π v_t

#A1_Y: Output Gap - Technology
function A1_Y(parameters::fixed_parameters)
    β=parameters.β
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π
    ρ_a=parameters.ρ_a
    σ=parameters.σ

    A1_coeff=-(-κ(parameters)* Ω(parameters)^2* (1-β*ϕ_π)* ψ_ya(parameters)* ρ_a* (ϕ_y+(1-ρ_a)*σ)
    -Ω(parameters)* ψ_ya(parameters)* (ϕ_y+(1-ρ_a)*σ)* (1-Ω(parameters)* ρ_a* (κ(parameters)+β*(ϕ_y+σ)) ) )/
    (κ(parameters)* Ω(parameters)^2* (1-β*ϕ_π)* ρ_a^2* σ - (1-Ω(parameters)*ρ_a*σ)* (1-Ω(parameters)* ρ_a*(κ(parameters) + β*(ϕ_y+σ))))
    
    return A1_coeff
end

#A1_π: Inflation - Technology
function A1_π(parameters::fixed_parameters)
    β=parameters.β
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π
    ρ_a=parameters.ρ_a
    σ=parameters.σ

    A1_coeff=-(κ(parameters)*Ω(parameters)*ϕ_y*ψ_ya(parameters) + κ(parameters)*Ω(parameters)*σ*ψ_ya(parameters)-κ(parameters)*Ω(parameters)*ρ_a*σ*ψ_ya(parameters))/
    (β*Ω(parameters)^2*σ^2*ρ_a^2 + β*κ(parameters)*Ω(parameters)^2*ϕ_π*σ*ρ_a^2 + β*Ω(parameters)^2*ϕ_y*σ*ρ_a^2-κ(parameters)*Ω(parameters)*ρ_a-β*Ω(parameters)*ϕ_y*ρ_a-β*Ω(parameters)*σ*ρ_a - Ω(parameters)*σ*ρ_a + 1)
    return A1_coeff
end

#A2_Y: Output Gap - Discount factor
function A2_Y(parameters::fixed_parameters)
    β=parameters.β
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π
    ρ_z=parameters.ρ_z
    σ=parameters.σ

    A2_coeff=-(κ(parameters)*Ω(parameters)^2*(1-β*ϕ_π)*(1-ρ_z)*ρ_z + Ω(parameters)*(1-ρ_z)*(1-Ω(parameters)*ρ_z*(κ(parameters)+β*(ϕ_y+σ))))/
    (κ(parameters)*Ω(parameters)^2*(1-β*ϕ_π)*ρ_z^2*σ-(1-Ω(parameters)*ρ_z*σ)*(1-Ω(parameters)*ρ_z*(κ(parameters)+β*(ϕ_y+σ))))
    
    return A2_coeff
end

#A2_π: Inflation - Discount factor
function A2_π(parameters::fixed_parameters)
    β=parameters.β
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π
    ρ_z=parameters.ρ_z
    σ=parameters.σ

    A2_coeff=-(-κ(parameters)*Ω(parameters)+κ(parameters)*Ω(parameters)*ρ_z)/
    (1-κ(parameters)*Ω(parameters)*ρ_z-β*Ω(parameters)*ϕ_y*ρ_z-Ω(parameters)*ρ_z*σ-β*Ω(parameters)*ρ_z*σ+β*κ(parameters)*Ω(parameters)^2*ϕ_π*ρ_z^2*σ+β*Ω(parameters)^2*ϕ_y*ρ_z^2*σ + β*Ω(parameters)^2*ρ_z^2*σ^2)    
    
    return A2_coeff
end


#A3_Y: Output Gap - Monetary Policy
function A3_Y(parameters::fixed_parameters)
    β=parameters.β
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π
    ρ_v=parameters.ρ_v
    σ=parameters.σ

    A3_coeff=-(-κ(parameters)*Ω(parameters)^2*(1-β*ϕ_π)*ρ_v-Ω(parameters)*(1-Ω(parameters)*ρ_v*(κ(parameters)+β*(ϕ_y+σ))))/
    (κ(parameters)*Ω(parameters)^2*(1-β*ϕ_π)*ρ_v^2*σ - (1-Ω(parameters)*ρ_v*σ)*(1-Ω(parameters)*ρ_v*(κ(parameters) + β*(ϕ_y+σ))))
    
    return A3_coeff
end

function A3_π(parameters::fixed_parameters)
    β=parameters.β
    ϕ_y=parameters.ϕ_y
    ϕ_π=parameters.ϕ_π
    ρ_v=parameters.ρ_v
    σ=parameters.σ

    A3_coeff=-(κ(parameters)*Ω(parameters))/
    (1-κ(parameters)*Ω(parameters)*ρ_v - β*Ω(parameters)*ϕ_y*ρ_v - Ω(parameters)*ρ_v*σ-β*Ω(parameters)*ρ_v*σ + β*κ(parameters)*Ω(parameters)^2*ϕ_π*ρ_v^2*σ  + β*Ω(parameters)^2*ϕ_y*ρ_v^2*σ + β*Ω(parameters)^2*ρ_v^2*σ^2)

    return A3_coeff
end

######################################################################################################
#Policy functions
#Policy function for inflation
#First dimension is a_t, second dimension is z_t, third dimension is v_t
function pol_output(grids::state_vars, coeff::log_lin_coeff)
    A1Y=coeff.A1Y
    A2Y=coeff.A2Y
    A3Y=coeff.A3Y

    a_grid=grids.a_t
    z_grid=grids.z_t
    v_grid=grids.v_t

    policy_function_output_gap=[A1Y*a_grid[i] + A2Y*z_grid[j] + A3Y*v_grid[k] for i in 1:length(a_grid), j in 1:length(z_grid), k in 1:length(v_grid)]

    return policy_function_output_gap
end

function pol_infl(grids::state_vars, coeff::log_lin_coeff)
    A1π=coeff.A1π
    A2π=coeff.A2π
    A3π=coeff.A3π

    a_grid=grids.a_t
    z_grid=grids.z_t
    v_grid=grids.v_t

    policy_function_inflation=[A1π*a_grid[i] + A2π*z_grid[j] + A3π*v_grid[k] for i in 1:length(a_grid), j in 1:length(z_grid), k in 1:length(v_grid)]

    return policy_function_inflation
end