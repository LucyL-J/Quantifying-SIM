using Distributions, SpecialFunctions, Optim, StatsBase

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

# Mutation rate estimation algorithms (if the optimisation fails, AIC=Inf is returned)
# Joint inference permissive+stressful condition
function estimate_mu_hom(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s; fitness_m=1.) # Estimating mutation rates under permissive/stressful conditons for the homogeneous-response model with optional differential fitness of mutants
    if fitness_m == "joint"
        mu_p, rho_p, AIC_p = estimate_mu_hom(mc_p, Nf_p, fitness_m="infer")
        mu_s, rho_s, AIC_s = estimate_mu_hom(mc_s, Nf_s, fitness_m="infer")
        log_likelihood_para_3(para) = -log_likelihood(mc_p, para[1], Nf_p, mc_s, Nf_s, para[2], para[3]) 
        res = Optim.optimize(log_likelihood_para_3, [mu_p, mu_s, rho_p])
        if Optim.converged(res) == true
            return[Optim.minimizer(res)[1],Optim.minimizer(res)[3],Optim.minimizer(res)[2],Optim.minimizer(res)[3], 6 + 2*Optim.minimum(res)]
        else
            return [0., -1., 0., -1., Inf]
        end
    else
        x_p = estimate_mu_hom(mc_p, Nf_p, fitness_m=fitness_m)                               # Estimation of mutation rate and optional differential fitness: permissive
        x_s = estimate_mu_hom(mc_s, Nf_s, fitness_m=fitness_m)                               # Estimation of mutation rate and optional differential fitness: stress
        return [x_p[1:end-1]; x_s[1:end-1]; x_p[end]+x_s[end]]                               # Returns inferred parameters plus AIC value (inferences permissive/stressful are independent)
    end
end

# Inference for a single condition
function estimate_mu_hom(mc::Vector{Int}, Nf; fitness_m=1.)                        # Estimating the mutation rate for a homogeneous population with optional differential fitness of mutants
    if fitness_m == "infer"                                                        # Differential fitness of mutants is inferred 
        log_likelihood_para_2(para) = -log_likelihood(mc, para[1], para[2], Nf)
        res = Optim.optimize(log_likelihood_para_2, [initial_mu(mc, Nf, 100), 1.]) # 2 inference parameters
        if Optim.converged(res) == true
            return [Optim.minimizer(res); 4 + 2*Optim.minimum(res)]                # Mutation rate, differential fitness of mutants, AIC
        else
            return [0., -1., Inf]
        end                                                   
    else                                                                           # Differential fitness of mutants is not inferred and set to 1 instead
        log_likelihood_para_1(para) = -log_likelihood(mc, para, fitness_m, Nf)
        res = Optim.optimize(log_likelihood_para_1, 0., maximum(mc)/Nf)            # 1 inference parameter
        if Optim.converged(res) == true
            return [Optim.minimizer(res), 2 + 2*Optim.minimum(res)]                # Mutation rate, AIC
        else
            return [0., Inf]                                                       
        end                                                                              
    end
end

# Log-likelihood functions
# Log-likelihood to observe a mutant count mc for a homogeneous population with differential fitness of mutants
function log_likelihood(mc::Vector{Int}, mutation_per_gen, fitness_m, Nf) 
    if mutation_per_gen<=0. || fitness_m<0. || Nf<=0.                     # Boundaries of the parameter regime
        return -Inf
    else
        p = P_mutant_count(maximum(mc), mutation_per_gen, Nf, fitness_m=fitness_m)
        mc[mc.>1000] .= 1000                                              # Any mutant count >1000 is considered as =1000 instead
        return sum(counts(mc, 0:maximum(mc)) .* log.(p))
    end
end
# Log-likelihood to observe mutant counts under permissive/stressful conditions with a joint differential fitness of mutants
function log_likelihood(mc_p::Vector{Int}, mutation_p_per_gen, Nf_p, mc_s::Vector{Int}, Nf_s, mutation_s_per_gen, fitness_m) 
    if mutation_p_per_gen<=0. || Nf_p<=0. || Nf_s<=0. || mutation_s_per_gen<=0. || fitness_m<0.                               # Boundaries of the parameter regime
        return -Inf
    else
        p_p = P_mutant_count(maximum(mc_p), mutation_p_per_gen, Nf_p, fitness_m=fitness_m)                                    # Mutant count distribution without stress
        p_s = P_mutant_count(maximum(mc_s), mutation_s_per_gen, Nf_s, fitness_m=fitness_m) 
        mc_p[mc_p.>1000] .= 1000                                                                                              # A mutant count >1000 is considered as =1000 instead
        mc_s[mc_s.>1000] .= 1000
        return sum(counts(mc_p, 0:maximum(mc_p)) .* log.(p_p)) + sum(counts(mc_s, 0:maximum(mc_s)) .* log.(p_s))              # The two observations are independent and their probabilities can be multiplied
    end
end

# Probability generating function method used to set initial value of the mutation rate based on
# Gillet-Markowska, A., Louvel, G., & Fischer, G. (2015). bz-rates: A web tool to estimate mutation rates from fluctuation analysis. G3: Genes, Genomes, Genetics, 5(11), 2323–2327. https://doi.org/10.1534/g3.115.019836
function empirical_pgf(z, x) # Empirical probability generating function calculated from a vector of observed data
    g = 0
    for i in x
        g += z^i
    end
    g /= length(x)
    return g
end
function initial_mu(z, mc::Vector{Int}, Nf)             # Estimate the mutation rate for a homogeneous population and given z        
    if z == 0.
        return log(empirical_pgf(z, mc)) / Nf
    else
        return z/((1-z)*log(1-z)) * log(empirical_pgf(z, mc)) / Nf
    end
end
function initial_mu(mc::Vector{Int}, Nf, z_values::Int) # Estimate the mutation rate for a homogeneous population by averaging over a number of z values
    mu = 0.
    for i = 0:z_values-1
        mu += initial_mu(i/z_values, mc, Nf)
    end
    return maximum([mu/z_values, 0.])
end