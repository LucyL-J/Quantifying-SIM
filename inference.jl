using SpecialFunctions, Optim, StatsBase

# The following recursive formulas for the pdf of the mutant count distribution are based on 
# Keller, P., & Antal, T. (2015). Mutant number distribution in an exponentially growing population. Journal of Statistical Mechanics: Theory and Experiment, 2015(1), P01011. https://doi.org/10.1088/1742-5468/2015/01/P01011

# Muntant count distribution for a homogeneous population with optional differential fitness of mutants
function P_mutant_count(K::Int, m; fit_m=1.)    # Mutations per generation = mutation rate [1/h] / growth rate [1/h]
    p = zeros(Float64, K+1)
    if fit_m == 0.
        for k = 0:K
            p[k+1] = m^k * exp(-m) / factorial(big(k))          # When the division rate is zero, the mutant count distribution is given by a Poisson distribution
        end
    else
        q = Q(K, m, fit_m)
        p[1] = exp(q[1])
        for k = 1:K
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
function Q(K::Int, m, fit_m)
    q = zeros(Float64, K+1)
    q[1] = -m
    if fit_m == 1.
        for k = 1:K
            q[k+1] = m / (k*(k+1))
        end
    else
        for k = 1:K
            q[k+1] = m/fit_m * factorial(big(k-1)) * gamma(1+1/fit_m) / gamma(1+1/fit_m+k)
        end
    end
    return q
end

# Mutant count distribution for a heterogeneous population with response-off and -on subpopulation. The relative division rate of on cells is an optional input parameter with default value of zero
# Mutation rates (for both response-off and -on cells) are given in mutations per division of response-off cells and scaled to be in units of mutations per generation
function P_mutant_count(K::Int, m, mu_het, f_on; rel_div_on=0.)
    if rel_div_on == 0.
        p_off = P_mutant_count(K, m)
        p_on = P_mutant_count(K, m * mu_het, fit_m=0.)
    else
        p_off = P_mutant_count(K, m * (1-f_on)/(1-f_on*(1-rel_div_on)))                                   
        p_on = P_mutant_count(K, m * mu_het * (1-f_on)/(1-f_on*(1-rel_div_on)), fit_m=rel_div_on/(1-f_on*(1-rel_div_on))) 
    end
    p = zeros(Float64, length(p_off))
    for k = 0:length(p_off)-1
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
# Estimating mutation rates under permissive/stressful conditons using the homogeneous-response model with optional differential fitness of mutants (jointly) inferred
function estimate_mu_hom(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s; fit_m=1.) 
    if fit_m == "joint"                                                                                  # Differential fitness of mutants is a joint inference parameter
        m_p, rho_p, AIC_p = estimate_init_hom(mc_p, fit_m="infer")                                       # Used as initial parameter in optimisation
        if AIC_p == Inf
            m_p, rho, AIC_p = estimate_init_hom(mc_p)
        end
        m_s, rho_s, AIC_s = estimate_init_hom(mc_s, fit_m="infer")                                       # Used as initial parameter in optimisation
        if AIC_s == Inf
            m_s, rho, AIC_s = estimate_init_hom(mc_s)
        end
        log_likelihood_para_3(para) = -log_likelihood(mc_p, para[1], mc_s, para[2], para[3])             # 3 inference parameters
        res = Optim.optimize(log_likelihood_para_3, [m_p, m_s, rho_s])                                   # Numbers of mutations permissive/stress, differential fitness
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            return[p[1]/Nf_p, p[3], p[2]/Nf_s, p[3], 6 + 2*Optim.minimum(res)]                           # Returns mutation rates permissive/stress, differential fitness plus AIC
        else
            return [0., -1., 0., -1., Inf]
        end
    else
        m_p, rho_p, AIC_p = estimate_init_hom(mc_p, fit_m=fit_m)              # Estimation of number of mutations and optional differential fitness: permissive
        m_s, rho_s, AIC_s = estimate_init_hom(mc_s, fit_m=fit_m)              # Estimation of number of mutations and optional differential fitness: stress
        return [m_p/Nf_p, rho_p, m_s/Nf_s, rho_s, AIC_p+AIC_s]                # Returns inferred parameters plus AIC value (inferences permissive/stressful are independent)
    end
end
# Estimating mutation rates of response-off/-on cells using the heterogeneous-response model with optional relative division rate of response-on cells inferred for unknown fraction of response-on subpopulation
function estimate_mu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s; rel_div_on=0.) 
    m = estimate_init_hom(mc_p)[1]*Nf_s/Nf_p                                                                          # Used as initial parameter in optimisation: number of mutations from response-off subpopulation
    mu_inc = estimate_init_hom(mc_s)[1] / m                                                                           # Used as initial parameter in optimisation: increase in population-average mutation rate for relative division rate of response-on cells = 1                                                                                           
    mu_het, div_init = estimate_init_het(mc_s, m, 0.)                                                                 # Used as initial parameter in optimisation: mutation-rate heterogeneity
    f_on = 1 - mu_inc/(mu_het+1)                                                                                      # Used as initial parameter in optimisation: final fraction of response-on subpopulation 
    if f_on <= 0.
        f_on = 1/mu_inc
    end
    f_on = - log(1-f_on) / log(Nf_s)                                                                                  # Used as initial parameter in optimisation: relative switching rate (equals fraction of response-on subpopulation for relative division rate of response-on cells = 0)
    log_likelihood_para_init(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], 0., para[3])            # Used as fixed parameter in optimisation: fraction of response-on subpopulation
    res = Optim.optimize(log_likelihood_para_init, [m, mu_het, f_on])
    if Optim.converged(res) == true
        p_init = Optim.minimizer(res)
        if rel_div_on == "infer"                                                                                      # Relative division rate of response-on cells is inferred 
            mu_het, rel_div_on = estimate_init_het(mc_s, p_init[1], p_init[3])                                        # Used as initial parameter in optimisation: mutation-rate heterogeneity, relative division rate of response-on subpopulation (initial value set to 0.5)
            f_on = p_init[3]/(1 - rel_div_on)                                                                         # Correct the fraction of response-on subpopulation by the relative division rate of response-on cells
            if f_on > 1. || f_on < 0.
                f_on = p_init[3] 
            end
            log_likelihood_para_3(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], para[3], f_on)     # 3 inference parameters (+1 from above)
            res = Optim.optimize(log_likelihood_para_3, [p_init[1], mu_het, rel_div_on])                              # Number of mutations stress, mutation-rate heterogeneity, relative division rate response-on cells 
            if Optim.converged(res) == true
                p = Optim.minimizer(res)
                return[p[1]/Nf_s, p[1]/Nf_s*p[2]*(1-f_on)/f_on, f_on, p[3], 8 + 2*Optim.minimum(res)]                 # Returns mutation rate response-off/-on cells, fraction of response-on subpopulation, relative division rate response-on cells, AIC
            else
                return [0., 0., 0., -1., Inf]
            end
        else
            mu_het, div_init = estimate_init_het(mc_s, p_init[1], p_init[3], rel_div_on=rel_div_on)                     # Used as initial parameter in optimisation: mutation-rate heterogeneity, used for correction of the fraction of response-on cells: relative division rate (initial value set to input parameter)
            f_on = p_init[3]/(1 - rel_div_on)                                                                           # Correct the fraction of response-on subpopulation by the inferred relative division rate of response-on cells 
            if f_on > 1. || f_on < 0.
                f_on = p_init[3] 
            end
            log_likelihood_para_2(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], rel_div_on, f_on)  # 2 inference parameters (+1 from above)
            res = Optim.optimize(log_likelihood_para_2, [p_init[1], p_init[2]])                                       # Number of mutations stress, mutation-rate heterogeneity
            if Optim.converged(res) == true
                p = Optim.minimizer(res)
                return[p[1]/Nf_s, p[1]/Nf_s*p[2]*(1-f_on)/f_on, f_on, rel_div_on, 6 + 2*Optim.minimum(res)]           # Returns mutation rate response-off/-on cells, fraction of response-on subpopulation, relative division rate response-on cells (input parameter), AIC
            else
                return [0., 0., 0., -1., Inf]
            end
        end
    else
        return [0., 0., 0., -1., Inf]
    end
end
# Estimating mutation rates of off/on cells using the new method with optional relative division rate of on cells inferred for known fraction of response-on subpopulation
function estimate_mu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, f_on; rel_div_on=0.)                  
    m = estimate_init_hom(mc_p)[1]*Nf_s/Nf_p                                                                     # Used as initial parameter in optimisation
    if rel_div_on == "infer"                                                                                     # Relative division rate of response-on cells is inferred
        mu_het, rel_div_on = estimate_init_het(mc_s, m, f_on)                                                    # Used as initial parameters in optimisation 
        log_likelihood_para_3(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], para[3], f_on)    # 3 inference parameters
        res = Optim.optimize(log_likelihood_para_3, [m, mu_het, rel_div_on])                                     # Number of mutations stress, mutation-rate heterogeneity, relative division rate response-on cells 
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            return [p[1]/Nf_s, p[2]*p[1]*(1-f_on)/(f_on*Nf_s), p[3], 6 + 2*Optim.minimum(res)]                   # Returns mutation rate response-off/-on cells, relative division rate response-on cells, AIC
        else
            return [0., 0., -1., Inf]                                                                               
        end
    else
        mu_het, div_init = estimate_init_het(mc_s, m, f_on, rel_div_on=rel_div_on)                               # Used as initial parameter in optimisation
        log_likelihood_para_2(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], rel_div_on, f_on) # 2 inference parameters
        res = Optim.optimize(log_likelihood_para_2, [m, mu_het])                                                 # Number of mutations stress, mutation-rate heterogeneity
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            return [p[1]/Nf_s, p[2]*p[1]*(1-f_on)/(f_on*Nf_s), rel_div_on, 4 + 2*Optim.minimum(res)]             # Returns mutation rate response-off/-on cells, relative division rate response-on cells (input parameter), AIC
        else
            return [0., 0., -1., Inf]                                                                                    
        end
    end
end


# The following inference functions are used to determine the initial parameters for the joint inference
function estimate_init_hom(mc::Vector{Int}; fit_m=1.)                          # Estimating the number of mutations for a homogeneous population with optional differential fitness of mutants
    if fit_m == "infer"                                                        # Differential fitness of mutants is inferred 
        log_likelihood_para_2(para) = -log_likelihood(mc, para[1], para[2])
        res = Optim.optimize(log_likelihood_para_2, [initial_m(mc, 100), 1.])  # 2 inference parameters
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            return [p[1], p[2], 4 + 2*Optim.minimum(res)]                      # Number of mutations, differential fitness of mutants, AIC
        else
            return [1., 1., Inf]
        end                                                   
    else                                                                       # Differential fitness of mutants is given (default set to 1)
        log_likelihood_para_1(para) = -log_likelihood(mc, para, fit_m)
        res = Optim.optimize(log_likelihood_para_1, 0., maximum(mc))           # 1 inference parameter
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            return [p, fit_m, 2 + 2*Optim.minimum(res)]                        # Number of mutations, differential fitness mutants (input parameter), AIC
        else
            return [1., 1., Inf]                                                       
        end                                                                              
    end
end
function estimate_init_het(mc::Vector{Int}, m, f_on; rel_div_on=0.)                                       # Estimation of the mutation-rate heterogeneity for given number of mutations (stress) and known fraction of response-on subpopulation
    if f_on == 0.
        log_likelihood_para_1(para) = -log_likelihood(mc, m, para, 0., 0.)                                # 1 inference parameter: mutation-rate heterogeneity
        res = Optim.optimize(log_likelihood_para_1, 0., maximum(mc))
        if Optim.converged(res) == true
            return [Optim.minimizer(res), rel_div_on]
        else
            return [maximum([initial_mu_het(mc, m, 100), f_on]), rel_div_on]
        end
    else
        log_likelihood_para_2(para) = -log_likelihood(mc, m, para[1], para[2], f_on)                      # 2 inference parameters: mutation-rate heterogeneity, relative division rate of response-on cells
        res = Optim.optimize(log_likelihood_para_2, [initial_mu_het(mc, m, 100), rel_div_on])             # Initial value for relative division rate of response-on cells is set input value
        if Optim.converged(res) == true
            return Optim.minimizer(res)                                                                   
        else
            return [maximum([initial_mu_het(mc, m, 100), f_on]), rel_div_on]
        end
    end
end

# Log-likelihood functions
# Log-likelihood to observe a mutant count mc for a homogeneous population with differential fitness of mutants
function log_likelihood(mc::Vector{Int}, m, fit_m) 
    if m<=0. || fit_m<0.                            # Boundaries of the parameter regime
        return -Inf
    else
        p = P_mutant_count(maximum(mc), m, fit_m=fit_m)
        return sum(counts(mc, 0:maximum(mc)) .* log.(p))
    end
end
# Joint log-likelihood to observe mutant counts under permissive/stressful conditions with a joint differential fitness of mutants
function log_likelihood(mc_p::Vector{Int}, m_p, mc_s::Vector{Int}, m_s, fit_m) 
    if m_p<=0. || m_s<=0. || fit_m<0.                                                              # Boundaries of the parameter regime
        return -Inf
    else
        p_p = P_mutant_count(maximum(mc_p), m_p, fit_m=fit_m)                                      # Mutant count distribution without stress
        p_s = P_mutant_count(maximum(mc_s), m_s, fit_m=fit_m) 
        return sum(counts(mc_p, 0:maximum(mc_p)) .* log.(p_p)) + sum(counts(mc_s, 0:maximum(mc_s)) .* log.(p_s))  # The two observations are independent and their probabilities can be multiplied
    end
end
# Log-likelihood to observe a mutant count mc for a heterogeneous population with relative division rate of response-on cells 
function log_likelihood(mc::Vector{Int}, m, mu_het, rel_div_on, f_on) 
    if m<=0. || mu_het<0. || rel_div_on<0. || f_on<0. || f_on>1.          # Boundaries of the parameter regime
        return -Inf
    else
        p = P_mutant_count(maximum(mc), m, mu_het, f_on, rel_div_on=rel_div_on)
        return sum(counts(mc, 0:maximum(mc)) .* log.(p))
    end
end
# Joint log-likelihood to observe the mutant counts (a) mc_p and (b) mc_s for 
# (a) homogeneous population without differential fitness of mutants
# (b) heterogeneous population with relative division rate of response-on cells 
function log_likelihood(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio, m, mu_het, rel_div_on, f_on)  
    if m<=0. || mu_het<0. || rel_div_on<0. || f_on<0. || f_on>1.                                                 # Boundaries of the parameter regime
        return -Inf
    else
        p_p = P_mutant_count(maximum(mc_p), m*N_ratio)                                                            # Mutant count distribution: permissive
        p_s = P_mutant_count(maximum(mc_s), m, mu_het, f_on, rel_div_on=rel_div_on)                               # Mutant count distribution: stress
        return sum(counts(mc_p, 0:maximum(mc_p)) .* log.(p_p)) + sum(counts(mc_s, 0:maximum(mc_s)) .* log.(p_s))  # The two observations are independent and their probabilities can be multiplied
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
function initial_m(z, mc::Vector{Int})                     # Estimate number of mutations for a homogeneous population and given z        
    if z == 0.
        return log(empirical_pgf(z, mc)) 
    else
        return z/((1-z)*log(1-z)) * log(empirical_pgf(z, mc))
    end
end
function initial_m(mc::Vector{Int}, z_values::Int)         # Estimate number of mutations for a homogeneous population by averaging over a number of z values
    m = 0.
    for i = 0:z_values-1
        m += initial_m(i/z_values, mc)
    end
    return maximum([m/z_values, 0.])
end
function initial_mu_het(z, mc::Vector{Int}, m)             # Estimate the mutation-rate heterogeneity under the heterogeneous-response model for given z   
    if z == 0.
        return -log(empirical_pgf(z, mc))/m - 1
    else
        return -log(empirical_pgf(z, mc))/(m*(1-z)) + log(1-z)/z
    end
end
function initial_mu_het(mc::Vector{Int}, m, z_values::Int) # Estimate the mutation-rate heterogeneity under the heterogeneous-response model by averaging over a number of z values 
    mu_het = 0.
    for i = 0:z_values-1
        mu_het += initial_mu_het(i/z_values, mc, m)
    end
    if mu_het == Inf || mu_het < 0.
        return 0.
    else
        return mu_het/z_values
    end
end