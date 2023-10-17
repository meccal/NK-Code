"""
LUCA MECCA
lmecca@london.edu
July 2023
"""

#These values are fixed
#From Judd (1998)

#This function computes the combinations of nodes and weights for the Chebyshev quadrature
#This function computes the Gauss-Hermite nodes for each collocation nodes of the Chebyshev polynomials
#Tensor product approach
function product_rule(constant_values::fixed_values, policy_f::policy_functions)::Tuple{Vector{Matrix{Float64}},Matrix{Float64}}
    state_vars=constant_values.state_variables
    quadrature_n=constant_values.quadrature_n #degree of approximation of the GH quadrature
    nodes=(constant_values.quadrature_nodes)[quadrature_n] #quadrature nodes
    weights=(constant_values.quadrature_weights)[quadrature_n] #quadrature weights
    state_n=constant_values.state_n #number of state variables
    collocation_n=constant_values.collocation_n #number of collocaiton points
    d_t=policy_f.d_t
    collocation_points=constant_values.collocation_points

    #create vector of nodes 
    nodes_vector=Vector{Vector{Vector{Float64}}}(undef, collocation_n)
    for i in 1:collocation_n
       nodes_vector[i]=Vector{Vector{Float64}}(undef, state_n)
        for j in 1:state_n
            current_state=state_vars[j]
            if current_state==d_lag
                #the endogenous state variable d_t requires special treatment becaus eit is not simply an AR(1) process
                #We know with certainty what will be the value of the endogenous state variable next period, there are no shocks
                nodes_vector[i][j]=state_to_Ch(current_state,ones(lastindex(nodes)).*d_t[i])
            else
                nodes_vector[i][j]=state_to_Ch(current_state, sqrt(2).*current_state.σ.*nodes.+ 
                (current_state.μ + current_state.ρ.*Ch_to_state(current_state,collocation_points[j,i])))
            end
        end
    end

    #vector of weights
    weights_vector=Vector{Vector{Float64}}(undef, state_n)
    for j in 1:state_n
        weights_vector[j]=weights
    end

    matrix_nodes=Vector{Matrix{Float64}}(undef, collocation_n)
    for i in 1:collocation_n
        tensor_nodes=[]
        tensor_zeros(nodes_vector[i], [], tensor_nodes)
        matrix_nodes[i] = hcat(tensor_nodes...)
    end

    tensor_weights=[]
    tensor_zeros(weights_vector, [], tensor_weights)
    matrix_weights = hcat(tensor_weights...)
   
    return matrix_nodes, matrix_weights
end


#This functions approximates the expectations term on the equilbrium condition for the first auxiliary variable 
#using the Gauss-Hermite quadrature
function GH_X1(constant_values::fixed_values, coefficients::coeffs, GH_elements::GH_nodes_weights)::Vector{Float64}
    
    GH_nodes_combination=GH_elements.nodes_combination
    GH_weights_combination=GH_elements.weights_combination
    ϵ=constant_values.ϵ
    α=constant_values.α
    collocation_n=constant_values.collocation_n #number of collocation points 
    state_n=constant_values.state_n #number of state variables
    coefficients_X1=coefficients.X1_t1
    coefficients_Π=coefficients.Π_t1

    #Final output is one expected value for each Chebyshev node
    expected_value=Vector{Float64}(undef, collocation_n)

    for i in 1:collocation_n
        #Redefine the constant values because of new collocation points
        constant_values_GH=fixed_values()
        constant_values_GH.state_variables, constant_values_GH.state_n=all_state_variables, lastindex(all_state_variables)
        constant_values_GH.collocation_points=GH_nodes_combination[i] #collocation points change at every iteration for quadrature
        constant_values_GH.collocation_n, constant_values_GH.combination_n=size(constant_values_GH.collocation_points,2), constant_values.combination_n
        constant_values_GH.Ch_bases_combinations=Ch_bases_function(constant_values_GH) #bases functions
        
        #evaluate the Chebyshev polynomial at the GH nodes
        Ch_polynomial=exp.(Ch_pol(constant_values_GH, coefficients_X1)).*
        exp.(Ch_pol(constant_values_GH, coefficients_Π)).^(ϵ/(1-α))
        
        
        expected_value[i]=π^(-0.5*state_n)*sum(transpose(prod(GH_weights_combination,dims=1)).*Ch_polynomial)
        
    end

    return expected_value

end


#This functions approximates the expectations term on the equilbrium condition for the second auxiliary variable 
#using the Gauss-Hermite quadrature
function GH_X2(constant_values::fixed_values, coefficients::coeffs, GH_elements::GH_nodes_weights, policy_f::policy_functions)::Vector{Float64}
    GH_nodes_combination=GH_elements.nodes_combination
    GH_weights_combination=GH_elements.weights_combination
    ϵ=constant_values.ϵ
    collocation_n=constant_values.collocation_n #number of collocation points 
    state_n=constant_values.state_n #number of state variables
    coefficients_X1=coefficients.X1_t1
    coefficients_Π=coefficients.Π_t1

    #Final output is one expected value for each Chebyshev node
    expected_value=Vector{Float64}(undef, collocation_n)

    for i in 1:collocation_n
        #Redefine the constant values because of new collocation points
        constant_values_GH=fixed_values()
        constant_values_GH.state_variables, constant_values_GH.state_n=all_state_variables, lastindex(all_state_variables)
        constant_values_GH.collocation_points=GH_nodes_combination[i] #collocation points change at every iteration for quadrature
        constant_values_GH.collocation_n, constant_values_GH.combination_n=size(constant_values_GH.collocation_points,2), constant_values.combination_n
        constant_values_GH.Ch_bases_combinations=Ch_bases_function(constant_values_GH) #bases functions
        
        #evaluate the Chebyshev polynomial at the GH nodes

        Ch_polynomial=x2_update(constant_values_GH,exp.(Ch_pol(constant_values_GH, coefficients_X1)), Π_star(constant_values, exp.(Ch_pol(constant_values_GH, coefficients_Π)))).*
        exp.(Ch_pol(constant_values_GH, coefficients_Π)).^(ϵ-1)     
        

        expected_value[i]=π^(-0.5*state_n)*sum(transpose(prod(GH_weights_combination,dims=1)).*Ch_polynomial)
    end

    return expected_value

end


#This functions approximates the expectations term on Euler equation using the Gauss-Hermite quadrature
function GH_EE(constant_values::fixed_values, coefficients::coeffs, GH_elements::GH_nodes_weights)::Vector{Float64}
    
    GH_nodes_combination=GH_elements.nodes_combination
    GH_weights_combination=GH_elements.weights_combination
    collocation_n=constant_values.collocation_n #number of collocation points 
    state_n=constant_values.state_n #number of state variables
    coefficients_C=coefficients.C_t1
    coefficients_Π=coefficients.Π_t1
    GH_nodes_combination=GH_elements.nodes_combination
    GH_weights_combination=GH_elements.weights_combination
    state_vars=constant_values.state_variables
    σ=constant_values.σ
   
    #Final output is one expected value for each Chebyshev node
    expected_value=Vector{Float64}(undef, collocation_n)

    for i in 1:collocation_n
        #Redefine the constant values because of new collocation points
        constant_values_GH=fixed_values()
        constant_values_GH.state_variables, constant_values_GH.state_n=all_state_variables, lastindex(all_state_variables)
        constant_values_GH.collocation_points=GH_nodes_combination[i] #collocation points change at every iteration for quadrature
        constant_values_GH.collocation_n, constant_values_GH.combination_n=size(constant_values_GH.collocation_points,2), constant_values.combination_n
        constant_values_GH.Ch_bases_combinations=Ch_bases_function(constant_values_GH) #bases functions
            
        #evaluate the Chebyshev polynomial at the GH nodes
        Ch_polynomial=exp.(Ch_to_state(state_vars[2],constant_values_GH.collocation_points[2,:]))./(exp.(Ch_pol(constant_values_GH, coefficients_C)).^σ.*exp.(Ch_pol(constant_values_GH, coefficients_Π)))
        
        expected_value[i]=π^(-0.5*state_n)*sum(transpose(prod(GH_weights_combination,dims=1)).*Ch_polynomial)
    end

    return expected_value

end