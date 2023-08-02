using Distributions, SpecialFunctions

# The following recursive formulas for the pdf of the mutant count distribution are based on 
# Keller, P., & Antal, T. (2015). Mutant number distribution in an exponentially growing population. Journal of Statistical Mechanics: Theory and Experiment, 2015(1), P01011. https://doi.org/10.1088/1742-5468/2015/01/P01011

# Muntant count distribution for a homogeneous population with optional differential fitness of mutants
function P_mutant_count(K::Int, mutation_per_gen, Nf; fitness_m=1.)    # Mutations generation = mutation rate [1/h] / growth rate [1/h]
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
