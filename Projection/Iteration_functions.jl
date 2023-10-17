"""
LUCA MECCA
lmecca@london.edu
July 2023
"""
#This file includes the iteration functions of the NK model

#Starting from an initial guess for the coefficients of the policy function Π_t, this function returns
#the coefficients and policy functions that make the two definitions of marginal cost equal
function iteration_Π_t(constant_values::fixed_values, policy_f::policy_functions, coefficients::coeffs)::Tuple{policy_functions, coeffs}
    #coefficients.Π_t=zeros(constant_values.collocation_n)

    function f!(F, ζ)
        
        coefficients.Π_t=ζ
      
        #a. Guess the policy function for Π_t
        policy_f.Π_t=exp.(Ch_pol(constant_values, coefficients.Π_t))

        #b. Compute Π_t^* from Π_t
        policy_f.Π_star_t=Π_star(constant_values,  policy_f.Π_t)

        #c. Compute d_t from its law of motion
        policy_f.d_t=law_motion_d(constant_values, policy_f)

        #d. Compute expectations
        GH_elements=GH_nodes_weights()
        GH_elements.nodes_combination,GH_elements.weights_combination=product_rule(constant_values, policy_f)
        policy_f.exp_X1=GH_X1(constant_values, coefficients, GH_elements) #expectation term of the equilibrium condition for the first auxiliary variable
        policy_f.exp_X2=GH_X2(constant_values, coefficients, GH_elements, policy_f) #expectation term of the equilibrium condition for the second auxiliary variable
        policy_f.exp_C=GH_EE(constant_values, coefficients, GH_elements) #expectations term of the Euler equation

        #e. Compute consumption from the Euler equation
        policy_f.C_t=C_EE(constant_values,policy_f)

        #f. Compute X_{2,t} from its intertemporal equation
        policy_f.X2_t=x2_int(constant_values,policy_f)

        #g. compute x_{1,t} from x_{2,t}
        policy_f.X1_t=x1_update(constant_values, policy_f.X2_t, policy_f.Π_star_t)
    
        #h. Compute Employment from the aggregate supply 
        policy_f.N_t=aggregate_supply(constant_values,policy_f)

        #i. Compute the real wage from the labour-leisure choice
        policy_f.real_wage=LL_choice(constant_values,policy_f)

        #compare the real marginal costs computed in two different ways
        policy_f.mc_1=mc_hat(constant_values, policy_f)
        policy_f.mc_2=mc(constant_values, policy_f)

        for i in 1:lastindex(policy_f.mc_1)
            F[i] = policy_f.mc_1[i]-policy_f.mc_2[i]
        end

    end

   
    sol=nlsolve(f!, coefficients.Π_t, method = :newton) #find the policy funciton Π_t that satisfied the conditions
    
    #redefine the coefficients and the policy functions
    coefficients.Π_t=sol.zero
    policy_f.Π_t=exp.(Ch_pol(constant_values, coefficients.Π_t)) #Policy function for Π_t
    policy_f.Π_star_t=Π_star(constant_values,  policy_f.Π_t)
    policy_f.d_t=law_motion_d(constant_values, policy_f)
    GH_elements=GH_nodes_weights()
    GH_elements.nodes_combination,GH_elements.weights_combination=product_rule(constant_values, policy_f)
    policy_f.exp_X1=GH_X1(constant_values, coefficients, GH_elements) #expectation term of the equilibrium condition for the first auxiliary variable
    policy_f.exp_X2=GH_X2(constant_values, coefficients, GH_elements, policy_f) #expectation term of the equilibrium condition for the second auxiliary variable
    policy_f.exp_C=GH_EE(constant_values, coefficients, GH_elements) #expectations term of the Euler equation
    policy_f.X2_t=x2_int(constant_values,policy_f)
    policy_f.X1_t=x1_update(constant_values, policy_f.X2_t, policy_f.Π_star_t)
    policy_f.C_t=C_EE(constant_values,policy_f)

    return policy_f, coefficients
    
end

#This function is the comprehensive solution algorithm of the model
#Algorithm steps
#1. Guess the policy functions for C, X1, Π and consider them as the true time t+1 functions
#2. a. Find the policy function for Π_t such that the two expressions for marginal cost coincide
#   b. Check if the functions C_t, X1_t and Π_t coincide with the time t+1 functions C_{t+1}, X1_{t+1}, Π_{t+1}.
#   If they coincide, stop the algorith, otherwise update the policy functions and coefficients and restart.

function solution_algorithm(constant_values::fixed_values, policy_f::policy_functions, coefficients::coeffs)::Tuple{policy_functions, coeffs}
    global difference=1
    global counter=0
    global error_list=ones(0)

    while difference>10^(-5) && counter<500
        global counter+=1
       
        print("Iteration number " * string(counter)*"\n")
        
        #Taking the guesses on the policy functions for C_{t+1}, X1_{t+1}, Π_{t+1} as true, this function finds the coefficients 
        #for Π_t such that the two definitions of marginal cost are equal
        @time global policy_f, coefficients=iteration_Π_t(constant_values, policy_f, coefficients)
       
        #compare the policy functions for C_{t+1}, X1_{t+1}, Π_{t+1} with their time-t counterparts
        #C_t, X1_t, Π_t
        global difference=sum(abs.(log.(policy_f.C_t1).-log.(policy_f.C_t))) + sum(abs.(log.(policy_f.X1_t1).-log.(policy_f.X1_t))) + sum(abs.(log.(policy_f.Π_t1).-log.(policy_f.Π_t)))
        append!(error_list, difference)
        print("The difference is " * string(difference) * string(counter)*"\n")

        #if the difference is greater than 10^(-6), then update the policy functions and the coefficients
        global policy_f.C_t1, policy_f.X1_t1, policy_f.Π_t1=copy(policy_f.C_t), copy(policy_f.X1_t), copy(policy_f.Π_t)
        bases=constant_values.Ch_bases_combinations
        global coefficients.C_t1=inv(bases' * bases)*bases'*log.(policy_f.C_t1)
        global coefficients.X1_t1=inv(bases' * bases)*bases'*log.(policy_f.X1_t1)
        global coefficients.Π_t1=inv(bases' * bases)*bases'*log.(policy_f.Π_t1)
    end

    return policy_f, coefficients

end