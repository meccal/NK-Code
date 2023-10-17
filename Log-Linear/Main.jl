"""
LUCA MECCA
lmecca@london.edu
Replicate the classical New Keynesian Model (Gali textbook, chapter 3)
LOG-LINEAR SOLUTION
July 2023
Written and tested in julia 1.8
"""

using Distributions, Plots

include("Struct.jl")
include("Discretize.jl")
include("Functions.jl")
include("Generic.jl")


#We have 3 state variables (in rlevant order for the code)
#1. Firm's technology
#2. Household's preference shock
#3. Monetary policy shock

#We use a simple interest rate rule for monetary policy (section 3.4.1 of Gali textbook)

################################################################
########################## PARAMETERS ##########################
################################################################
#Parameters are from Chapter 3 of the Gali textbook
#Household
const σ=1 #CRRA coefficient
const ϕ=5 #Frisch elasticity of labour supply
const β=0.99 #ssubjective discount rate
const ϵ=9 #elasticity of substitution between consumption goods

const ρ_z=0.5 #persistence of preference shock
const σ_z=1 #volatility of preference shock
const n_z=15 #number of grid points for the preference shock


#Firm
const α=0.25 #curvature of labour
const θ=0.75 #probability of not resetting prices

const ρ_a=0.9 #persistence of TFP shock
const σ_a=1 #volatility of TFP shock
const n_a=15 #number of grid points for the TFP shock


#Taylor rule
const ϕ_π=1.5 #MP coefficient on inflation
const ϕ_y=0.5/4 #MP coefficient on output

const ρ_v=0.5 #persistence of MP shock
const σ_v=1 #volatility of MP shock
const n_v=15 #number of grid points for the MP shock

#Monedy Demand
const η=4 #semielasticity of money demand

#IRF
const qrts=15 #number of quarters of the IRF
const ϵ_v=0.25 #size of the MP shock
const ϵ_z=-0.5 #size of the preference shock
const ϵ_a=1 #size of the tech shock

#Define the parameters data structure
parameters=fixed_parameters(σ,ϕ,β,ϵ,ρ_z,σ_z,n_z,α,θ,ρ_a,σ_a, n_a,ϕ_π,ϕ_y,ρ_v,σ_v,n_v,η,qrts,ϵ_v,ϵ_z,ϵ_a)



################################################################################
########################## IMPULSE RESPONSE FUNCTIONS ##########################
################################################################################
#Monetary policy shocks (section 3.4.1.1)
MP_shocks=[parameters.σ_v*parameters.ϵ_v*parameters.ρ_v^i for i in 0:parameters.qrts] #compute the MP shocks
state=state_vars(zeros(parameters.qrts+1), zeros(parameters.qrts+1), MP_shocks) #define the state variables
#recall that for a monetary policy shock, z_t and a_t are assumed to be 0 at all t

IRF_output_gap_MP=IRF_MP_Y(parameters, state) #Effect on output (gap)
#It coincides with output because the natural level of outut is unchanged after a MP shock
IRF_inflation_MP=IRF_MP_π(parameters, state, "Yes") #Effect on inflation (annualized)
IRF_real_rate_MP=IRF_MP_r(parameters, state, "Yes") #Effect on the real interest rate (annualized)
IRF_nominal_rate_MP=IRF_MP_i(parameters, state, "Yes") #Effect on the nominal interest rate (annualized)
IRF_employment_MP=empl(parameters, state, IRF_output_gap_MP) #Effect on employment
IRF_real_wage_MP=real_wage(parameters, IRF_output_gap_MP, IRF_employment_MP) #Effect on the real wage
IRF_price_level_MP=cumsum(IRF_MP_π(parameters, state)) #Effect on the price level
IRF_money_supply_MP=IRF_MP_Ms(parameters, state, IRF_price_level_MP)


################################################################################
############################### POLICY FUNCTIONS ###############################
################################################################################

###############################################################
#1. Define grids for the state variables
#We will need the grids to plot the policy functions. We do not need the grid for the simulation
#We use the simple Tauchen methodology
a_grid, _=Tau(ρ_a, σ_a, n_a) #technology shock
z_grid, _=Tau(ρ_z, σ_z, n_z) #preference shock
v_grid, _=Tau(ρ_v, σ_v, n_v) #MP shock

grids=state_vars(a_grid, z_grid, v_grid)

##############################################################################################################################
#2. Compute log-linear coefficients
coeff=log_lin_coeff(A1_Y(parameters),A2_Y(parameters),A3_Y(parameters), A1_π(parameters),A2_π(parameters),A3_π(parameters))

##############################################################################################################################
#3. Define policy functions for inflation and output gap
policy_function_output_gap=pol_output(grids, coeff) #output gap
policy_function_inflation=pol_infl(grids, coeff) #inflation