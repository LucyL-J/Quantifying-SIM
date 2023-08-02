using Distributions, SpecialFunctions

# The following recursive formulas for the pdf of the mutant count distribution are based on 
# Keller, P., & Antal, T. (2015). Mutant number distribution in an exponentially growing population. Journal of Statistical Mechanics: Theory and Experiment, 2015(1), P01011. https://doi.org/10.1088/1742-5468/2015/01/P01011

# Muntant count distribution for a homogeneous population with optional differential fitness of mutants
function P_mutant_count(K::Int, mutation_per_gen, Nf; fitness_m=1.)    # Mutations per generation = mutation rate [1/h] / growth rate [1/h]
    K_bound = minimum([K, 1000])                                       # Introducing an arbritraty cut-off at 1000 colonies; any mutant count >1000 is considered as =1000 instead
    p = zeros(Float64, K_bound+1)
    if fitness_m == 0.
        for k = 0:K_bound
            p[k+1] = pdf(Poisson(Nf*mutation_per_gen), k)              # When the division rate is zero, the mutant count distribution is given by a Poisson distribution
        end
    else
        q = Q(K_bound, mutation_per_gen, fitness_m, Nf)
        p[1] = exp(q[1])
        for k = 1:K_bound
            S = 0.
            for i = 0:k-1
                S += (k-i) * q[k-i+1] * p[i+1]
            end
            p[k+1] = S/k
        end
    end
    return p
end

# Recursive helper functions 
function Q(K::Int, mutation_per_gen, fitness_m, Nf)
    q = zeros(Float64, K+1)
    q[1] = -Nf*mutation_per_gen
    if fitness_m == 1.
        for k = 1:K
            q[k+1] = Nf*mutation_per_gen / (k*(k+1))
        end
    else
        for k = 1:K
            q[k+1] = Nf*mutation_per_gen/fitness_m * factorial(big(k-1)) * gamma(1+1/fitness_m) / gamma(1+1/fitness_m+k)
        end
    end
    return q
end

# Mutant count distribution for a heterogeneous population with response-off and -on subpopulation. The relative division rate of on cells is an optional input parameter with default value of zero
# Mutation rates (for both response-off and -on cells) are given in mutations per division of response-off cells and scaled to be in units of mutations per generation
function P_mutant_count(K::Int, mutation_off_per_div, Nf, mutation_on_per_div, f_on; rel_division_on=0.)
    p_off = P_mutant_count(K, mutation_off_per_div/(1-f_on*(1-rel_division_on)), (1-f_on)*Nf)                                   
    p_on = P_mutant_count(K, mutation_on_per_div/(1-f_on*(1-rel_division_on)), f_on*Nf, fitness_m=rel_division_on/(1-f_on*(1-rel_division_on))) 
    K_bound = length(p_off)-1
    p = zeros(Float64, K_bound+1)
    for k = 0:K_bound
        pk = 0
        for j = 0:k
            pk += p_off[j+1] * p_on[k-j+1] # pdf of the total mutant count (response-off + -on)
        end
        p[k+1] = pk
    end
    return p
end