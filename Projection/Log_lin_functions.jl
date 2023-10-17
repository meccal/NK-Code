"""
LUCA MECCA
lmecca@london.edu
July 2023
"""

#This file contains the functions used to compute the log-linear policy functions and impulse response functions

#OPTIMAL MARKUP
function μ(constant_values::fixed_values)::Float64
    ϵ=constant_values.ϵ

    opt_markup=log(ϵ/(ϵ-1))
    return opt_markup
end

#RESPONSE OF INFLATION TO THE OUTPUT GAP
function κ(constant_values::fixed_values)::Float64
    θ=constant_values.θ
    α=constant_values.α
    ϵ=constant_values.ϵ
    β=constant_values.β
    σ=constant_values.σ
    ϕ=constant_values.ϕ

    kappa=(1-θ)*(1-β*θ)/θ*((1-α)/(1-α+α*ϵ))*(σ+(ϕ+α)/(1-α))

    return kappa

end

function Ω(constant_values::fixed_values)::Float64
    σ=constant_values.σ
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π  

    Omega=1/(σ+ϕ_y+κ(constant_values)*ϕ_π)
    return Omega
end


function ψ_ya(constant_values::fixed_values)::Float64
    σ=constant_values.σ
    ϕ=constant_values.ϕ
    α=constant_values.α  

    psi_ya=(1+ϕ)/(σ*(1-α)+ϕ+α)
    return psi_ya
end


function ψ_y(constant_values::fixed_values)::Float64
    σ=constant_values.σ
    ϕ=constant_values.ϕ
    α=constant_values.α  

    psi_y=-((1-α)*(μ(constant_values)-log(1-α)))/(σ*(1-α)+ϕ+α)
    return psi_y
end


#EMPLOYMENT
#y_t is output
#a_t is the technology shock
function empl(constant_values::fixed_values, a_t::Vector{Float64}, y_t::Vector{Float64})::Vector{Float64}
    α=constant_values.α

    n_t=(y_t-a_t)./(1-α)
    return n_t
end


#REAL WAGE
#y_t is output
#a_t is the technology shock
function real_wage(constant_values::fixed_values, y_t::Vector{Float64}, n_t::Vector{Float64})::Vector{Float64}
    σ=constant_values.σ
    ϕ=constant_values.ϕ

    w_t_real=σ.*y_t+ϕ.*n_t
    return w_t_real
end


#Effect of monetary policy shock on output (gap)
#annual:input "Yes" if you want the IRF to be annualized, "No" otherwise
function IRF_MP_Y(constant_values::fixed_values, shock::Vector{Float64})::Vector{Float64}
    β=constant_values.β
    ρ_v=constant_values.state_variables[3].ρ   
    σ=constant_values.σ
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π

    MP_shocks=shock

    Y_MP=(-(1-β*ρ_v)/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(constant_values)*(ϕ_π-ρ_v))).*MP_shocks

    return Y_MP
end


#Effect of monetary policy shock on inflation
#annual:input "Yes" if you want the IRF to be annualized, "No" otherwise
function IRF_MP_π(constant_values::fixed_values, shock::Vector{Float64})::Vector{Float64}
    β=constant_values.β
    ρ_v=constant_values.state_variables[3].ρ   
    σ=constant_values.σ
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π

    MP_shocks=shock
  
    π_MP=(-κ(constant_values)/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(constant_values)*(ϕ_π-ρ_v))).*MP_shocks.*4 #annualized

    return π_MP
end


#Effect of monetary policy shock on the real interest rate
#annual:input "Yes" if you want the IRF to be annualized, "No" otherwise
function IRF_MP_r(constant_values::fixed_values, shock::Vector{Float64})::Vector{Float64}
    β=constant_values.β
    ρ_v=constant_values.state_variables[3].ρ   
    σ=constant_values.σ
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π

    MP_shocks=shock
   
    r_MP=(σ*(1-ρ_v)*(1-β*ρ_v)/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(constant_values)*(ϕ_π-ρ_v))).*MP_shocks.*4 #annualized

    return r_MP
end


#Effect of monetary policy shock on the nominal interest rate
#annual: input "Yes" if you want the IRF to be annualized, "No" otherwise
function IRF_MP_i(constant_values::fixed_values, shock::Vector{Float64})::Vector{Float64}
    β=constant_values.β
    ρ_v=constant_values.state_variables[3].ρ   
    σ=constant_values.σ
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π

    MP_shocks=shock
   
    i_MP=((σ*(1-ρ_v)*(1-β*ρ_v)-ρ_v*κ(constant_values))/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(constant_values)*(ϕ_π-ρ_v))).*MP_shocks.*4 #annualized
    

    return i_MP
end


#Effect of monetary policy shock on money supply
#annual: input "Yes" if you want the IRF to be annualized, "No" otherwise
#p_t is the price level
function IRF_MP_Ms(constant_values::fixed_values, shock::Vector{Float64}, p_t::Vector{Float64})::Vector{Float64}
    β=constant_values.β
    ρ_v=constant_values.state_variables[3].ρ   
    σ=constant_values.σ
    ϕ_y=constant_values.ϕ_y
    η=constant_values.η
    ϕ_π=constant_values.ϕ_π

    MP_shocks=shock
    lag_price=lag(p_t,1)
    lag_price[1]=0
   
    Ms_MP=lag_price.-(((1-β*ρ_v)*(1+η*σ*(1-ρ_v))+(1-η*ρ_v)*κ(constant_values))/((1-β*ρ_v)*(σ*(1-ρ_v)+ϕ_y) +κ(constant_values)*(ϕ_π-ρ_v))).*MP_shocks

    return Ms_MP
end


#This function computes the log-linear IRFs of different variables to a MP shock
function MP_IRF(constant_values::fixed_values)
    σ_v=constant_values.state_variables[3].σ
    ρ_v=constant_values.state_variables[3].ρ
    qrts=constant_values.qrts 
    ϵ_v=constant_values.ϵ_v 
    
    #Monetary policy shocks (section 3.4.1.1)
    MP_shocks=[σ_v*ϵ_v*ρ_v^i for i in 0:qrts] #compute the MP shocks
    a_t=zeros(constant_values.qrts+1)

    #recall that for a monetary policy shock, z_t and a_t are assumed to be 0 at all t
    IRF_output_gap_MP=IRF_MP_Y(constant_values, MP_shocks) #Effect on output (gap)
    #It coincides with output because the natural level of outut is unchanged after a MP shock
    IRF_inflation_MP=IRF_MP_π(constant_values,MP_shocks) #Effect on inflation (annualized)
    IRF_real_rate_MP=IRF_MP_r(constant_values, MP_shocks) #Effect on the real interest rate (annualized)
    IRF_nominal_rate_MP=IRF_MP_i(constant_values, MP_shocks) #Effect on the nominal interest rate (annualized)
    IRF_employment_MP=empl(constant_values, a_t, IRF_output_gap_MP) #Effect on employment
    IRF_real_wage_MP=real_wage(constant_values, IRF_output_gap_MP, IRF_employment_MP) #Effect on the real wage
    IRF_price_level_MP=cumsum(IRF_MP_π(constant_values, MP_shocks)) #Effect on the price level
    IRF_money_supply_MP=IRF_MP_Ms(constant_values, MP_shocks, IRF_price_level_MP)

    return IRF_output_gap_MP, IRF_inflation_MP, IRF_real_rate_MP, IRF_nominal_rate_MP, IRF_employment_MP, 
    IRF_real_wage_MP, IRF_price_level_MP, IRF_money_supply_MP, MP_shocks

end


######################################################################################################
#The following functions compute the coefficients of the following guesses for the policy functions:
#y_t = ψ_y + A1_y a_t + A2_y z_t + A3_y v_t
#π_t = A1_π a_t + A2_π z_t + A3_π v_t

#A1_Y: Output - Technology
function A1_Y(constant_values::fixed_values)
    β=constant_values.β
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π
    ρ_a=constant_values.state_variables[1].ρ
    σ=constant_values.σ

    A1_coeff=-(-Ω(constant_values)*(1-β*ϕ_π)*ρ_a*(κ(constant_values)*Ω(constant_values)*ψ_ya(constant_values)*ρ_a*σ+κ(constant_values)*Ω(constant_values)*ψ_ya(constant_values)*(ϕ_y+(1-ρ_a)*σ))-
    (-ψ_ya(constant_values)+Ω(constant_values)*ψ_ya(constant_values)*ρ_a*σ+Ω(constant_values)*ψ_ya(constant_values)*(ϕ_y+(1-ρ_a)*σ))*
    (1-Ω(constant_values)*ρ_a*(κ(constant_values)+β*(ϕ_y+σ))))/
    (κ(constant_values)* Ω(constant_values)^2* (1-β*ϕ_π)* ρ_a^2* σ - (1-Ω(constant_values)*ρ_a*σ)* (1-Ω(constant_values)* ρ_a*(κ(constant_values) + β*(ϕ_y+σ))))
    
    return A1_coeff
end

#A1_π: Inflation - Technology
function A1_π(constant_values::fixed_values)
    β=constant_values.β
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π
    ρ_a=constant_values.state_variables[1].ρ
    σ=constant_values.σ

    A1_coeff=(Ω(constant_values)*(-κ(constant_values)*ϕ_y*ψ_ya(constant_values)-κ(constant_values)*ψ_ya(constant_values)*σ+κ(constant_values)*ψ_ya(constant_values)*ρ_a*σ))/
    (1-κ(constant_values)*Ω(constant_values)*ρ_a-β*Ω(constant_values)*ϕ_y*ρ_a-Ω(constant_values)*ρ_a*σ-β*Ω(constant_values)*ρ_a*σ+β*κ(constant_values)*Ω(constant_values)^2*ϕ_π*ρ_a^2*σ+
    β*Ω(constant_values)^2*ϕ_y*ρ_a^2*σ+β*Ω(constant_values)^2*ρ_a^2*σ^2)
       
    return A1_coeff
end

#A2_Y: Output - Discount factor
function A2_Y(constant_values::fixed_values)
    β=constant_values.β
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π
    ρ_z=constant_values.state_variables[2].ρ
    σ=constant_values.σ

    A2_coeff=-(κ(constant_values)*Ω(constant_values)^2*(1-β*ϕ_π)*(1-ρ_z)*ρ_z + Ω(constant_values)*(1-ρ_z)*(1-Ω(constant_values)*ρ_z*(κ(constant_values)+β*(ϕ_y+σ))))/
    (κ(constant_values)*Ω(constant_values)^2*(1-β*ϕ_π)*ρ_z^2*σ-(1-Ω(constant_values)*ρ_z*σ)*(1-Ω(constant_values)*ρ_z*(κ(constant_values)+β*(ϕ_y+σ))))
    
    return A2_coeff
end

#A2_π: Inflation - Discount factor
function A2_π(constant_values::fixed_values)
    β=constant_values.β
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π
    ρ_z=constant_values.state_variables[2].ρ
    σ=constant_values.σ

    A2_coeff=-(-κ(constant_values)*Ω(constant_values)+κ(constant_values)*Ω(constant_values)*ρ_z)/
    (1-κ(constant_values)*Ω(constant_values)*ρ_z-β*Ω(constant_values)*ϕ_y*ρ_z-Ω(constant_values)*ρ_z*σ-β*Ω(constant_values)*ρ_z*σ+β*κ(constant_values)*Ω(constant_values)^2*ϕ_π*ρ_z^2*σ+β*Ω(constant_values)^2*ϕ_y*ρ_z^2*σ + β*Ω(constant_values)^2*ρ_z^2*σ^2)    
    
    return A2_coeff
end


#A3_Y: Output - Monetary Policy
function A3_Y(constant_values::fixed_values)
    β=constant_values.β
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π
    ρ_v=constant_values.state_variables[3].ρ
    σ=constant_values.σ

    A3_coeff=-(-κ(constant_values)*Ω(constant_values)^2*(1-β*ϕ_π)*ρ_v-Ω(constant_values)*(1-Ω(constant_values)*ρ_v*(κ(constant_values)+β*(ϕ_y+σ))))/
    (κ(constant_values)*Ω(constant_values)^2*(1-β*ϕ_π)*ρ_v^2*σ - (1-Ω(constant_values)*ρ_v*σ)*(1-Ω(constant_values)*ρ_v*(κ(constant_values) + β*(ϕ_y+σ))))
    
    return A3_coeff
end

function A3_π(constant_values::fixed_values)
    β=constant_values.β
    ϕ_y=constant_values.ϕ_y
    ϕ_π=constant_values.ϕ_π
    ρ_v=constant_values.state_variables[3].ρ
    σ=constant_values.σ

    A3_coeff=-(κ(constant_values)*Ω(constant_values))/
    (1-κ(constant_values)*Ω(constant_values)*ρ_v - β*Ω(constant_values)*ϕ_y*ρ_v - Ω(constant_values)*ρ_v*σ-β*Ω(constant_values)*ρ_v*σ + β*κ(constant_values)*Ω(constant_values)^2*ϕ_π*ρ_v^2*σ  + β*Ω(constant_values)^2*ϕ_y*ρ_v^2*σ + β*Ω(constant_values)^2*ρ_v^2*σ^2)

    return A3_coeff
end

#This function computes the coefficients of the log-linear policy functions
function log_lin_coefficients(constant_values::fixed_values)
    A1_coeff_Y=A1_Y(constant_values)
    A2_coeff_Y=A2_Y(constant_values)
    A3_coeff_Y=A3_Y(constant_values)
    A1_coeff_π=A1_π(constant_values)
    A2_coeff_π=A2_π(constant_values)
    A3_coeff_π=A3_π(constant_values)

    return A1_coeff_Y, A2_coeff_Y, A3_coeff_Y, A1_coeff_π, A2_coeff_π, A3_coeff_π
end 

#This function evaluates the log-linear policy functions at a defined set of points
#evaluation_points is the set of points the policy funciton needs to be evaluated at (not necessarily the collocaiton points in this case)
#policy allows to choose whether to compute the policy function for output or inflation
#input 'output' if you want to compute the policy function for output and 'inflation' otherwise
function policy_log_lin(eval_points::Matrix{Float64}, policy::String="output")::Vector{Float64}
    if !(policy in ["output", "inflation"])
        error("Invalid policy type. Allowed values are 'output' and 'inflation'.")
    end
    
    #compute the coefficients
    A1_coeff_Y, A2_coeff_Y, A3_coeff_Y, A1_coeff_π, A2_coeff_π, A3_coeff_π=log_lin_coefficients(constant_values)

    if policy=="output"
        A0=ψ_y(constant_values)
        coefficients= [A1_coeff_Y, A2_coeff_Y, A3_coeff_Y]
    else
        A0=0
        coefficients= [A1_coeff_π, A2_coeff_π, A3_coeff_π]
    end

    #Take the evaluation points, eliminate the last row 
    #(corresponding to the price dispersion, which is zero in the log-linear solution)
    #and convert the points from the Chebyshev space to the dominion of the corresponding state variable
    eval_points_log_lin=eval_points[1:end-1,:]

    for i in 1:size(eval_points_log_lin,1)
        state=constant_values.state_variables[i]
        eval_points_log_lin[i,:]=Ch_to_state(state, eval_points_log_lin[i,:])
    end

    policy_f_log_lin=exp.(A0.+ sum(coefficients.*eval_points_log_lin, dims=1))
 
    return vec(policy_f_log_lin)

end
