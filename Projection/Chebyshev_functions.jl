"""
LUCA MECCA
lmecca@london.edu
July 2023
"""

#This file contains functions related to Chebyshev polynomials

#This functions defines the hypercube state variables
#Takes as input a vector of state variables
#returns the boundaries computed as 3*unconditional standard deviation

function hypercube(state_vars::Vector{state_var})::Vector{state_var}

    for i in 1:lastindex(state_vars) #repeat the same process for each state variable
        
        if state_vars[i].lb==nothing && state_vars[i].ub==nothing #check if the boundaries are not predefined
            state_vars[i].lb, state_vars[i].ub=-3*state_vars[i].σ/sqrt(1-state_vars[i].ρ), 3*state_vars[i].σ/sqrt(1-state_vars[i].ρ)
        end

    end
    return state_vars
end


#This function transforms arguments defined over the Chebyshev domain [-1, 1] into outputs belonging to the domain of the state variable, defined by a min and a max
#x can be either a single number or a vector 
function Ch_to_state(state::state_var, x::Union{Float64, Int64, Vector{Float64}, Vector{Int64}})
        x_transformed=state.lb .+(0.5.*(state.ub-state.lb)).*(1 .+ x)
    return x_transformed
end


#This function transforms arguments defined the domain of the state variable into outputs belonging to the Chebyshev domain [-1, 1]
#x can be either a single number or a vector 
function state_to_Ch(state::state_var, x::Union{Float64, Int64, Vector{Float64}, Vector{Int64}})
       x_transformed=(2 .*x.-(state.lb+state.ub))./(state.ub-state.lb)
    return x_transformed
end


#This functions finds the zeros of the Chebyshev polynomials
#n is the degree of approximation of the Chebyshev polynomials
function Ch_zero(n::Int64=3)::Vector{Float64}
    Ch_zeros=Vector{Float64}(undef,n+1)
    for j in 1:n+1
        Ch_zeros[j]=cos(π*((2*j-1)/(2*(n+1))))
    end
    return Ch_zeros
end


#This recursive function takes as inputs a vector of vectors containing the Chebyshev zeros and finds all possible tensor combinations
function tensor_zeros(zeros_matrix::Vector{Vector{Float64}}, current_combination, tensor_combinations)
    if isempty(zeros_matrix)
        push!(tensor_combinations, copy(current_combination))
    else
        for element in zeros_matrix[1]
            push!(current_combination, element)
            tensor_zeros(zeros_matrix[2:end], current_combination, tensor_combinations)
            pop!(current_combination)
        end
    end
end


#This recursive function takes as inputs a vector of vectors containing the Chebyshev zeros and finds all possible combinations that do not exceed a predefined order of approximation
#p: maximum order of approximation
function complete_zeros(zeros_matrix::Vector{Vector{Float64}}, current_combination, complete_combinations, index_sum::Int64=0, p::Int64=3)
    if isempty(zeros_matrix)
        if index_sum <= p
            push!(complete_combinations, copy(current_combination))
        end
    else
        for (index, element) in enumerate(zeros_matrix[1])
            push!(current_combination, element)
            complete_zeros(zeros_matrix[2:end], current_combination, complete_combinations, index_sum + index-1,p)
            pop!(current_combination)
        end
    end
end


#This function takes as input a vector of vector of Chebyshev bases for single state variables and computes the tensor product
#This function takes as input a vector of vector of Chebyshev bases for single state variables and computes the tensor product
function tensor_product(Ch_base::Vector{Vector{Float64}})::Vector{Float64}
    if length(Ch_base) == 0
        return [1]
    end
    
    tensor_comb = []
    
    for element in Ch_base[1]
        sub_products = tensor_product(Ch_base[2:end])
        for sub_product in sub_products
            push!(tensor_comb, element * sub_product)
        end
    end
    
    return tensor_comb
end


#This function takes as input a vector of vector of Chebyshev bases for single state variables and computes the complete product
function complete_product(Ch_base::Vector{Vector{Float64}}, p::Int64=3)::Tuple{Vector{Float64}, Vector{Int64}}
    if length(Ch_base) == 0
        return [1], [0]
    end
    
    complete_comb = []
    index_comb=[]
    
    for i in 1:lastindex(Ch_base[1])
        element=Ch_base[1][i]
        sub_products, indices = complete_product(Ch_base[2:end],p)
        for j in 1:lastindex(sub_products)
            sub_product=sub_products[j]
            index=indices[j]
            push!(complete_comb, element * sub_product)
            push!(index_comb, i-1+index)
        end
    end

    #drop if the order of approximation exceeds p
    complete_comb=complete_comb[index_comb.<=p]
    index_comb=index_comb[index_comb.<=p]
    
    return complete_comb, index_comb
end


#This function computes all the combinations of Chebyshev zeros
#Takes as input a vector of state variables
#product allows to choose two options: choose "tensor" to have the n-fold tensor product as basis functions
#choose "complete" to have the complete set of polynomials of degree p.
function Ch_zero_combine(constant_values::fixed_values)::Matrix{Float64}
    product=constant_values.product
    state_vars=constant_values.state_variables
    if !(product in ["tensor", "complete"])
        error("Invalid product type. Allowed values are 'tensor' and 'complete'.")
    end

    #create a vector of vectors containing the Chebyshev zeros for each state variable
    Ch_zeros=Vector{Vector{Float64}}(undef,lastindex(state_vars))
    for i in 1:lastindex(state_vars)
        Ch_zeros[i]=Ch_zero(state_vars[i].n)
    end

    if product =="tensor"
        #create all the possible tensor combinations
        tensor_combinations=[]
        tensor_zeros(Ch_zeros, [], tensor_combinations)
        matrix = hcat(tensor_combinations...)

    else
        #If instead of the tensor product, we want the complete set, we need to eliminate some combinations in the columns of matrix
        p=maximum([state_vars[i].n for i in 1:lastindex(state_vars)]) #unique order of approximation for the complete set
       
        #creat all combinations that do not exceed order of approximation pairs
        complete_combinations=[]
        complete_zeros(Ch_zeros, [], complete_combinations,0,p)
        matrix = hcat(complete_combinations...)
    end

    return matrix
end


#This functions computes the Chebyshev bases (which are invariant to the different iterations) at all collocation points
#3. product allows to choose two options: choose "tensor" to have the n-fold tensor product as basis functions
#choose "complete" to have the complete set of polynomials of degree p.
function Ch_bases_function(constant_values::fixed_values)::Matrix{Float64}

    collocation_points=constant_values.collocation_points
    product=constant_values.product
    state_vars=constant_values.state_variables
    combination_n=constant_values.combination_n #number of combinations
    collocation_n=constant_values.collocation_n #number of collocation points (in quadrature can differ from the number of combinations)
    state_n=constant_values.state_n #number of state variables

    if !(product in ["tensor", "complete"])
        error("Invalid product type. Allowed values are 'tensor' and 'complete'.")
    end

    #This structure is quite complex
    #First layer: Each element of the most external vector stands for a combination of collocation points
    #For each collocation point, we have a vector of vectors containing the bases for each state variable, evaluated at that collocation point
    Ch_bases=Vector{Any}(undef, collocation_n)
    for i in 1:collocation_n
        Ch_bases[i]=Vector{Vector{Float64}}(undef,state_n)
    end

    for i in 1:lastindex(state_vars)
        current_state_variable=state_vars[i] #select the relevant state variable
        degree_approx=current_state_variable.n #and its degree of approximation

        Ch_base=Matrix{Float64}(undef,degree_approx+1, collocation_n)
        #First two bases
        Ch_base[1,:].=1
        Ch_base[2,:].=collocation_points[i,:]
        
        if degree_approx>1
        #Complete the bases using the recursive structure of Chebyshev polynomials
            for j in 3:degree_approx+1
                Ch_base[j,:]=2 .*collocation_points[i,:].*Ch_base[j-1,:].-Ch_base[j-2,:]
            end
        else
        end

        #Now place each base in the vector that corresponds to the relevant node
        for col in 1:collocation_n
            Ch_bases[col][i]=Ch_base[:,col]
        end
    end

    if product=="tensor"
        #Find the tensor combinations of the Chebyshev bases
        Ch_bases_combinations=Matrix{Float64}(undef,collocation_n,combination_n)
        for i in 1:collocation_n
            Ch_bases_combinations[i,:]=tensor_product(Ch_bases[i])
        end
    
    else
        #creat all combinations that do not exceed order of approximation pairs
        p=maximum([state_vars[i].n for i in 1:lastindex(state_vars)]) #unique order of approximation for the complete set
        Ch_bases_combinations=Matrix{Float64}(undef,collocation_n,combination_n)
        for i in 1:collocation_n
            Ch_bases_combinations[i,:]=complete_product(Ch_bases[i],p)[1]
        end
    end

    return Ch_bases_combinations
end


#This function approximates a polcy function using Chebyshev polynomials
#Inputs are:
function Ch_pol(constant_values::fixed_values,  coefficients::Vector{Float64})::Union{Float64, Vector{Float64}}
    Ch_bases_combinations=constant_values.Ch_bases_combinations

    policy_function=sum(transpose(coefficients).*Ch_bases_combinations, dims=2)
    return vec(policy_function)
end


#This function evaluates the policy functions at a defined set of points using Chebyshev polynomials
#eval_points is the set of points the policy funciton needs to be evaluated at (not necessarily the collocaiton points in this case)
#policy allows to choose whether to compute the policy function for output or inflation
#input 'output' if you want to compute the policy function for output and 'inflation' otherwise
function policy_Ch(constant_values::fixed_values, final_coefficients::coeffs, eval_points::Matrix{Float64}, policy::String="output")::Vector{Float64}
    
    if !(policy in ["output", "inflation"])
        error("Invalid policy type. Allowed values are 'output' and 'inflation'.")
    end

    #redefine constant values to incorporate the new evaluation points
    constant_values_new=fixed_values()
    constant_values_new.state_variables, constant_values_new.state_n=all_state_variables, lastindex(all_state_variables)
    constant_values_new.collocation_points=eval_points
    constant_values_new.collocation_n, constant_values_new.combination_n=size(constant_values_new.collocation_points,2), constant_values.combination_n
    constant_values_new.Ch_bases_combinations=Ch_bases_function(constant_values_new) #bases functions

    #Choose the coefficients depending on which policy function you want to compute
    if policy=="output"
        coefficients=final_coefficients.C_t1
    else
        coefficients=final_coefficients.Π_t1
    end

    #compute the policy function
    policy_f_Ch=exp.(Ch_pol(constant_values_new, coefficients))
    
    return policy_f_Ch

end


