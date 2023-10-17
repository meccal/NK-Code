"""
LUCA MECCA
lmecca@london.edu
Replicate the classical New Keynesian Model (Gali textbook 2015, chapter 3)
PROJECTION METHODS (COLLOCATION) USING CHEBYSHEV POLYNOMIALS
July 2023
Written and tested in julia 1.8
"""

using Distributions, Plots, Random, Combinatorics, NLsolve, MAT, Kronecker

include("Struct.jl")
include("Chebyshev_functions.jl")
include("Functions.jl")
include("Quadrature.jl")
include("Iteration_functions.jl")
include("Generic.jl")
include("Log_lin_functions.jl")


#We have 4 state variables (in relevant order for the code)
#1. Firm's technology (a_t)
#2. Household's preference shock (z_t)
#3. Monetary policy shock (ν_t)
#4. Price dispersion (D_{t-1})

#We use a simple interest rate rule for monetary policy (section 3.4.1 of Gali textbook)

################################################################
############################ SETUP #############################
################################################################
#All other parameters should be modified in the cosntructor fixed_values in the Struct.jl file
#Technology
const ρ_a=0.9
const σ_a=0.00025
const n_a=2 #order of approximation for technology

#Preference
const ρ_z=0.5
const σ_z=0.00025
const n_z=2 #order of approximation for preferences

#Monetary Policy
const ρ_v=0.5
const σ_v=0.00025
const n_v=2 #order of approximation for monetary policy

#Price dispersion
#Bounds for the endogenous state variables are from Fernandez-Villaverde, Gordon, Guerron, Rubio-Ramirez (2015)
const d_min=1
const d_max=1.005
const n_d=1 #order of approximation for price dispersion

#Define the state variables of the model
#constant_values collects all the values that are invariante to the iterations of the algorithm
constant_values=fixed_values()
a_t=state_var(ρ=ρ_a,σ=σ_a,n=n_a) #Technology shock
z_t=state_var(ρ=ρ_z,σ=σ_z,n=n_z) #Preference shock
v_t=state_var(ρ=ρ_v,σ=σ_v,n=n_v) #Monetary policy shock
d_lag=state_var(lb=d_min, ub=d_max, n=n_d) #Price dispersion
#include all the state variables in this vector
all_state_variables=[a_t, z_t, v_t, d_lag]
#Define the hypercube for the state variables and add the extremes to the state_variables
hypercube(all_state_variables)
constant_values.state_variables, constant_values.state_n=all_state_variables, lastindex(all_state_variables)
constant_values.collocation_points=Ch_zero_combine(constant_values) #collocation points
constant_values.collocation_n, constant_values.combination_n=size(constant_values.collocation_points,2), size(constant_values.collocation_points,2)
constant_values.Ch_bases_combinations=Ch_bases_function(constant_values) #bases functions
constant_values.steady_state=compute_SS(constant_values) #steady state of output

########################################################################################################################################################

################################################################
########################### SOLUTION ###########################
################################################################
#Policy functions guessed
#log(\hat{C}_t)=\hat{f^1}(a_t, z_t, v_t, d_{t-1})
#log(\hat{X1}_t)=\hat{f^2}(a_t, z_t, v_t, d_{t-1})
#log(\hat{Π}_t)=\hat{f^3}(a_t, z_t, v_t, d_{t-1})

#Initial guesses
#Initial coefficients are chosen such that the initial value of the policy functions are close to their steady state
coefficients=coeffs([log(constant_values.steady_state.Y);zeros(constant_values.collocation_n-1)], [log(constant_values.steady_state.X1);zeros(constant_values.collocation_n-1)], zeros(constant_values.collocation_n), zeros(constant_values.collocation_n))
#compute the log-linear coefficients
#log_coeffs=log_lin_coefficients(constant_values)
#coefficients=coeffs([ψ_y(constant_values), 0, Ch_to_state(all_state_variables[3], log_coeffs[3]), Ch_to_state(all_state_variables[2], log_coeffs[2]), Ch_to_state(all_state_variables[1], log_coeffs[1])], [log(constant_values.steady_state.X1);zeros(constant_values.collocation_n-1)], [0, 0, Ch_to_state(all_state_variables[3], log_coeffs[6]), Ch_to_state(all_state_variables[2],log_coeffs[5]), Ch_to_state(all_state_variables[2],log_coeffs[4])], [0, 0, Ch_to_state(all_state_variables[3], log_coeffs[6]), Ch_to_state(all_state_variables[2],log_coeffs[5]), Ch_to_state(all_state_variables[2],log_coeffs[4])] )
policy_f=policy_functions(exp.(Ch_pol(constant_values, coefficients.C_t1)), exp.(Ch_pol(constant_values, coefficients.X1_t1)), exp.(Ch_pol(constant_values, coefficients.Π_t1)) )

#Find solution
final_policy_f, final_coefficients=solution_algorithm(constant_values, policy_f, coefficients)

#Compare the log-linear and projection results
include("Comparison.jl")

