using Distributions, SpecialFunctions, Optim, StatsBase

# The following recursive formulas for the pdf of the mutant count distribution are based on 
# Keller, P., & Antal, T. (2015). Mutant number distribution in an exponentially growing population. Journal of Statistical Mechanics: Theory and Experiment, 2015(1), P01011. https://doi.org/10.1088/1742-5468/2015/01/P01011

# Muntant count distribution for a homogeneous population with optional differential fitness of mutants
function P_mutant_count(K::Int, m; fit_m=1.)    # Mutations per generation = mutation rate [1/h] / growth rate [1/h]
    p = zeros(Float64, K+1)
    if fit_m == 0.
        for k = 0:K
            p[k+1] = pdf(Poisson(m), k)          # When the division rate is zero, the mutant count distribution is given by a Poisson distribution
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
    p_off = P_mutant_count(K, m * (1-f_on)/(1-f_on*(1-rel_div_on)))                                   
    p_on = P_mutant_count(K, m * mu_het * (1-f_on)/(1-f_on*(1-rel_div_on)), fit_m=rel_div_on/(1-f_on*(1-rel_div_on))) 
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

# Mutation rate estimation and model selection between heterogeneous/homogeneous response
# Input parameters
# mc_p: Mutant counts for permissive condition
# Nf_p: Final population size for permissive condition
# mc_s: Mutant counts for stressful condition
# Nf_s: Final population size for stressful condition
# Optional
# fitm_p: Mutant fitness under permissive condition
# fitm_s: Mutant fitness under stressful condition
# To constrain the mutant fitness to be equal under permissive and stressful conditions, set joint=true
# Output: selected model and inferred parameters
function estimate_mu(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s; fitm_p=false, fitm_s=false, joint=false)
    mu_off, mu_on, f_on, AIC_het = estimate_mu_het(mc_p, Nf_p, mc_s, Nf_s)
    if joint == true
        mu_p_set, mu_s_set, AIC_set = estimate_mu_hom(mc_p, Nf_p, mc_s, Nf_s)
        mu_p_joint, mu_s_joint, rho_joint, AIC_joint = estimate_mu_hom(mc_p, Nf_p, mc_s, Nf_s, fit_m="joint")
        AIC_infer = Inf
    elseif fitm_p == false && fitm_s == false
        mu_p_set, mu_s_set, AIC_set = estimate_mu_hom(mc_p, Nf_p, mc_s, Nf_s)
        mu_p_joint, mu_s_joint, rho_joint, AIC_joint = estimate_mu_hom(mc_p, Nf_p, mc_s, Nf_s, fit_m="joint")
        mu_p_infer, rho_p_infer, mu_s_infer, rho_s_infer, AIC_infer = estimate_mu_hom(mc_p, Nf_p, mc_s, Nf_s, fit_m="infer")
    elseif fitm_s == false 
        AIC_set = Inf
        m_p_joint, AIC_p = estimate_init_hom(mc_p, fit_m=fitm_p)
        mu_p_joint = m_p_joint/Nf_p
        mu_p_infer = mu_p_joint
        m_s_joint, AIC_s_joint = estimate_init_hom(mc_s, fit_m=fitm_p)
        mu_s_joint = m_s_joint/Nf_s
        AIC_joint = AIC_p + AIC_s_joint
        m_s_infer, rho_s_infer, AIC_s_infer = estimate_init_hom(mc_s, fit_m="infer")
        mu_s_infer = m_s_infer/Nf_s
        AIC_infer = AIC_p + AIC_s_infer
        rho_p_infer = fitm_p
    elseif fitm_p == false 
        AIC_set = Inf
        m_s_joint, AIC_s = estimate_init_hom(mc_s, fit_m=fitm_s)
        mu_s_joint = m_s_joint/Nf_s
        mu_s_infer = mu_s_joint
        m_p_joint, AIC_p_joint = estimate_init_hom(mc_p, fit_m=fitm_s)
        mu_p_joint = m_p_joint/Nf_p
        AIC_joint = AIC_p_joint + AIC_s
        m_p_infer, rho_p, AIC_p_infer = estimate_init_hom(mc_p, fit_m="infer")
        mu_p_infer = m_p_infer/Nf_p
        AIC_infer = AIC_p_infer + AIC_s
        rho_s_infer = fitm_s  
    else
        m_p_set, AIC_p_set = estimate_init_hom(mc_p, fit_m=fitm_p)
        mu_p_set = m_p_set/Nf_p
        m_s_set, AIC_s_set = estimate_init_hom(mc_s, fit_m=fitm_s)
        mu_s_set = m_s_set/Nf_s
        AIC_set = AIC_p_set + AIC_s_set
        AIC_joint = Inf
        AIC_infer = Inf
    end
    AIC_hom = minimum([AIC_set, AIC_joint, AIC_infer])
    if AIC_het - AIC_hom < -2
        println("Heterogeneous-response model with non-dividing response-on cells is selected")
        println("Mutation rate response-off cells = ", mu_off)
        println("Mutation rate response-on cells = ", mu_on)
        println("Fraction response-on subpopulation = ", f_on)
        println("(Beta) Relative switching rate = ", (log(1-f_on)-log(f_on))/log(Nf_s))
        println("Relative mutation-rate increase = ", mu_on/mu_off)
        println("Increase in population mean mutation rate = ", (mu_off*(1-f_on) + mu_on*f_on)/mu_off)
    else
        if AIC_het - AIC_hom > 2
            println("Homogeneous-response model is selected")
        else
            println("No preferred model")
            println("")
            println("Heterogeneous-response model inferred parameters:")
            println("Mutation rate response-off cells = ", mu_off)
            println("Mutation rate response-on cells = ", mu_on)
            println("Fraction response-on subpopulation = ", f_on)
            println("(Beta) Relative switching rate = ", (log(1-f_on)-log(f_on))/log(Nf_s))
            println("Relative mutation-rate increase = ", mu_on/mu_off)
            println("Increase in population mean mutation rate = ", (mu_off*(1-f_on) + mu_on*f_on)/mu_off)
            println("")
            println("Homogeneous-response model inferred parameters:")
        end
        m = argmin([AIC_set, AIC_joint, AIC_infer])
        if m == 1
            println("(Mutant fitness set to 1 or the input value)")
            println("Mutation rate permissive condition = ", mu_p_set)
            println("Mutation rate stressful condition = ", mu_s_set)
            println("Increase in population mean mutation rate = ", mu_s_set/mu_p_set)
        elseif m == 2
            println("(Mutant fitness inferred, constrained to be equal under permissive and stressful conditions)")
            println("Mutation rate permissive condition = ", mu_p_joint)
            println("Mutation rate stressful condition = ", mu_s_joint)
            println("Mutant fitness = ", rho_joint)
            println("Increase in population mean mutation rate = ", mu_s_joint/mu_p_joint)
        elseif m == 3
            println("(Mutant fitness inferred)")
            println("Mutation rate permissive condition = ", mu_p_infer)
            println("Mutation rate stressful condition = ", mu_s_infer)
            println("Mutant fitness permissive condition = ", rho_p_infer)
            println("Mutant fitness stressful condition = ", rho_s_infer)
            println("Increase in population mean mutation rate = ", mu_s_infer/mu_p_infer)
        end
    end
end

# Mutation rate estimation algorithms (if the optimisation fails, AIC=Inf is returned)
# Joint inference permissive+stressful condition
# Estimating mutation rates under permissive/stressful conditons using the homogeneous-response model with optional differential fitness of mutants (jointly) inferred
function estimate_mu_hom(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s; fit_m=1.) 
    if fit_m == "joint"                                                                                  # Differential fitness of mutants is a joint inference parameter
        m_p, rho_p, AIC_p = estimate_init_hom(mc_p, fit_m="infer")                                       # Used as initial parameter in optimisation
        m_s, rho_s, AIC_s = estimate_init_hom(mc_s, fit_m="infer")                                       # Used as initial parameter in optimisation
        log_likelihood_para_3(para) = -log_likelihood(mc_p, para[1], mc_s, para[2], para[3])             # 3 inference parameters
        res = Optim.optimize(log_likelihood_para_3, [m_p, m_s, rho_s])                                   # Numbers of mutations permissive/stress, differential fitness
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            return[p[1]/Nf_p, p[2]/Nf_s, p[3], 6 + 2*Optim.minimum(res)]                                 # Returns mutation rates permissive/stress, differential fitness plus AIC
        else
            return [0., 0., -1., Inf]
        end
    else
        x_p = estimate_init_hom(mc_p, fit_m=fit_m)                            # Estimation of number of mutations and optional differential fitness: permissive
        x_p[1] /= Nf_p                                                        # Mutation rate = number of mutations / final population size
        x_s = estimate_init_hom(mc_s, fit_m=fit_m)                            # Estimation of number of mutations and optional differential fitness: stress
        x_s[1] /= Nf_s                                                        # Mutation rate = number of mutations / final population size
        return [x_p[1:end-1]; x_s[1:end-1]; x_p[end]+x_s[end]]                # Returns inferred parameters plus AIC value (inferences permissive/stressful are independent)
    end
end
# Estimating mutation rates of response-off/-on cells using the heterogeneous-response model with optional relative division rate of response-on cells for unknown fraction of response-on subpopulation
function estimate_mu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s) 
    m = estimate_init_hom(mc_p)[1]*Nf_s/Nf_p                                                                          # Used as initial parameter in optimisation
    mu_inc = estimate_init_hom(mc_s)[1] / m                                                                           # Used as initial parameters in optimisation                                                                                           
    mu_het = maximum([initial_mu_het(mc_s, m, 100), m/Nf_p])                                                          # Used as initial parameter in optimisation
    f_on = maximum([1 - mu_inc/(mu_het+1), 0.])                                                                       # Used as initial parameter in optimisation
    switching = log((1-f_on)/f_on) / log(Nf_s)
    log_likelihood_para_3(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], 0., para[3])               # 3 inference parameters
    res = Optim.optimize(log_likelihood_para_3, [m, mu_het, switching])                                               # Number of mutations stress, mutation-rate heterogeneity, relative switching rate          
    if Optim.converged(res) == true
        p = Optim.minimizer(res)
        return [p[1]/Nf_s, p[2]*p[1]*(Nf_s^(p[3]-1)), 1/(Nf_s^p[3]), 6 + 2*Optim.minimum(res)]                        # Returns mutation rate response-off/-on cells, fraction of response-on subpopulation, AIC
    else
        println(mu_het, f_on, switching)
        return [0., 0., 0., Inf]                                                                                                        
    end
end
function estimate_mu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, f_on; rel_div_on=0.)                  # Estimating mutation rates of off/on cells using the new method with optional relative division rate of on cells for known fraction of response-on subpopulation
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
        mu_het = maximum([initial_mu_het(mc_s, m, 100), m/Nf_p])                                                 # Used as initial parameter in optimisation
        log_likelihood_para_2(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], rel_div_on, f_on) # 2 inference parameters
        res = Optim.optimize(log_likelihood_para_2, [m, mu_het])                                                 # Number of mutations stress, mutation-rate heterogeneity
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            return [p[1]/Nf_s, p[2]*p[1]*(1-f_on)/(f_on*Nf_s), 4 + 2*Optim.minimum(res)]                         # Returns mutation rate response-off/-on cells, AIC
        else
            return [0., 0., Inf]                                                                                    
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
            return [p, 2 + 2*Optim.minimum(res)]                               # Number of mutations, AIC
        else
            return [1., Inf]                                                       
        end                                                                              
    end
end
function estimate_init_het(mc::Vector{Int}, m, f_on)                                                  # Estimation of the mutation-rate heterogeneity for given number of mutations (stress) and known fraction of response-on subpopulation
    log_likelihood_para_2(para) = -log_likelihood(mc, m, para[1], para[2], f_on)                      # 2 inference parameters: mutation-rate heterogeneity, relative division rate of response-on cells
    res = Optim.optimize(log_likelihood_para_2, [initial_mu_het(mc, m, 100), 0.5])                    # Initial value for relative division rate of response-on cells is set to 0.5
    if Optim.converged(res) == true
        return Optim.minimizer(res)                                                                   
    else
        return [1., 0.]
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
    return maximum([mu_het/z_values, 0.])
end