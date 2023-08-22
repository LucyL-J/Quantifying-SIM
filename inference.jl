using Distributions, SpecialFunctions, Optim, StatsBase

# The following recursive formulas for the pdf of the mutant count distribution are based on 
# Keller, P., & Antal, T. (2015). Mutant number distribution in an exponentially growing population. Journal of Statistical Mechanics: Theory and Experiment, 2015(1), P01011. https://doi.org/10.1088/1742-5468/2015/01/P01011

# Muntant count distribution for a homogeneous population with optional differential fitness of mutants
function P_mutant_count(K::Int, mu_per_gen, Nf; fit_m=1.)    # Mutations per generation = mutation rate [1/h] / growth rate [1/h]
    p = zeros(Float64, K+1)
    if fit_m == 0.
        for k = 0:K
            p[k+1] = pdf(Poisson(Nf*mu_per_gen), k)          # When the division rate is zero, the mutant count distribution is given by a Poisson distribution
        end
    else
        q = Q(K, mu_per_gen, fit_m, Nf)
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
function Q(K::Int, mu_per_gen, fit_m, Nf)
    q = zeros(Float64, K+1)
    q[1] = -Nf*mu_per_gen
    if fit_m == 1.
        for k = 1:K
            q[k+1] = Nf*mu_per_gen / (k*(k+1))
        end
    else
        for k = 1:K
            q[k+1] = Nf*mu_per_gen/fit_m * factorial(big(k-1)) * gamma(1+1/fit_m) / gamma(1+1/fit_m+k)
        end
    end
    return q
end

# Mutant count distribution for a heterogeneous population with response-off and -on subpopulation. The relative division rate of on cells is an optional input parameter with default value of zero
# Mutation rates (for both response-off and -on cells) are given in mutations per division of response-off cells and scaled to be in units of mutations per generation
function P_mutant_count(K::Int, mu_off_per_div, Nf, mu_on_per_div, f_on; rel_div_on=0.)
    p_off = P_mutant_count(K, mu_off_per_div/(1-f_on*(1-rel_div_on)), (1-f_on)*Nf)                                   
    p_on = P_mutant_count(K, mu_on_per_div/(1-f_on*(1-rel_div_on)), f_on*Nf, fit_m=rel_div_on/(1-f_on*(1-rel_div_on))) 
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
        mu_p_joint, AIC_p = estimate_mu_hom(mc_p, Nf_p, fit_m=fitm_p)
        mu_p_infer = mu_p_joint
        mu_s_joint, AIC_s_joint = estimate_mu_hom(mc_s, Nf_s, fit_m=fitm_p)
        AIC_joint = AIC_p + AIC_s_joint
        mu_s_infer, rho_s_infer, AIC_s_infer = estimate_mu_hom(mc_s, Nf_s, fit_m="infer")
        AIC_infer = AIC_p + AIC_s_infer
        rho_p_infer = fitm_p
    elseif fitm_p == false 
        AIC_set = Inf
        mu_s_joint, AIC_s = estimate_mu_hom(mc_s, Nf_s, fit_m=fitm_s)
        mu_s_infer = mu_s_joint
        mu_p_joint, AIC_p_joint = estimate_mu_hom(mc_p, Nf_p, fit_m=fitm_s)
        AIC_joint = AIC_p_joint + AIC_s
        mu_p_infer, rho_p, AIC_p_infer = estimate_mu_hom(mc_p, Nf_p, fit_m="infer")
        AIC_infer = AIC_p_infer + AIC_s
        rho_s_infer = fitm_s  
    else
        mu_p_set, AIC_p_set = estimate_mu_hom(mc_p, Nf_p, fit_m=fitm_p)
        mu_s_set, AIC_s_set = estimate_mu_hom(mc_s, Nf_s, fit_m=fitm_s)
        AIC_set = AIC_p_set + AIC_s_set
        AIC_joint = Inf
        AIC_infer = Inf
    end
    AIC_hom = minimum([AIC_set, AIC_joint, AIC_infer])
    if AIC_het - AIC_hom < -2
        println("Heterogeneous-response model is selected")
        println("Mutation rate response-off = ", mu_off)
        println("Mutation rate response-on = ", mu_on)
        println("Fraction of response-on subpopulation = ", f_on)
        println("Relative mutation-rate increase = ", mu_on/mu_off)
        println("Increase in population mean mutation rate = ", (mu_off*(1-f_on) + mu_on*f_on)/mu_off)
    else
        if AIC_het - AIC_hom > 2
            println("Homogeneous-response model is selected")
        else
            println("No preferred model")
            println("Heterogeneous-response model inferred parameters:")
            println("Mutation rate response-off = ", mu_off)
            println("Mutation rate response-on = ", mu_on)
            println("Fraction of response-on subpopulation = ", f_on)
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
        mu_p, rho_p, AIC_p = estimate_mu_hom(mc_p, Nf_p, fit_m="infer")                                  # Used as initial parameter in optimisation
        mu_s, rho_s, AIC_s = estimate_mu_hom(mc_s, Nf_s, fit_m="infer")                                  # Used as initial parameter in optimisation
        log_likelihood_para_3(para) = -log_likelihood(mc_p, para[1], Nf_p, mc_s, Nf_s, para[2], para[3]) # 3 inference parameters
        res = Optim.optimize(log_likelihood_para_3, [mu_p, mu_s, rho_p])                                 # Mutation rate permissive, mutation rate stress, differential fitness
        if Optim.converged(res) == true
            return[Optim.minimizer(res); 6 + 2*Optim.minimum(res)] 
        else
            return [0., 0., -1., Inf]
        end
    else
        x_p = estimate_mu_hom(mc_p, Nf_p, fit_m=fit_m)                               # Estimation of mutation rate and optional differential fitness: permissive
        x_s = estimate_mu_hom(mc_s, Nf_s, fit_m=fit_m)                               # Estimation of mutation rate and optional differential fitness: stress
        return [x_p[1:end-1]; x_s[1:end-1]; x_p[end]+x_s[end]]                       # Returns inferred parameters plus AIC value (inferences permissive/stressful are independent)
    end
end
# Estimating mutation rates of response-off/-on cells using the heterogeneous-response model with optional relative division rate of response-on cells for unknown fraction of response-on subpopulation
function estimate_mu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s; rel_div_on=0.) 
    mu_off = estimate_mu_hom(mc_p, Nf_p)[1]                                                                               # Used as initial parameter in optimisation
    if rel_div_on == "infer"                                                                                              # Relative division rate of response-on cells is inferred  
        mu_on, rel_div_on, f_on = estimate_mu_het(mc_s, mu_off, Nf_s, infer_div_on=true)                                  # Used as initial parameters in optimisation
        log_likelihood_para_4(para) = -log_likelihood(mc_p, para[1], Nf_p, mc_s, Nf_s, para[2], para[3], para[4])
        res = Optim.optimize(log_likelihood_para_4, [mu_off, mu_on, rel_div_on, f_on])                                    # 4 inference parameters
        if Optim.converged(res) == true
            return [Optim.minimizer(res); 8 + 2*Optim.minimum(res)]                                                       # Mutation rate response-off/-on cells, relative division rate response-on cells, fraction of response-on subpopulation, AIC
        else
            return [0., 0., -1., 0., Inf]                                                                                                                                                  
        end
    else                                                                                                                  # Relative division rate of response-on cells is given (default set to 0)
        mu_on, f_on = estimate_mu_het(mc_s, mu_off, Nf_s)                                                                # Used as initial parameters in optimisation
        log_likelihood_para_3(para) = -log_likelihood(mc_p, para[1], Nf_p, mc_s, Nf_s, para[2], rel_div_on, para[3])
        res = Optim.optimize(log_likelihood_para_3, [mu_off, mu_on, f_on])                                                # 3 inference parameters
        if Optim.converged(res) == true
            return [Optim.minimizer(res); 6 + 2*Optim.minimum(res)]                                                       # Mutation rate response-off/-on cells, fraction of response-on subpopulation, AIC
        else
            return [0., 0., 0., Inf]                                                                                                        
        end
    end
end
function estimate_mu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, f_on; rel_div_on=0.)   # Estimating mutation rates of off/on cells using the new method with optional relative division rate of on cells for known fraction of response-on subpopulation
    mu_off = estimate_mu_hom(mc_p, Nf_p)[1]                                                       # Used as initial parameter in optimisation
    if rel_div_on == "infer"                                                                      # Relative division rate of response-on cells is inferred
        mu_on, rel_div_on = estimate_mu_het(mc_s, mu_off, Nf_s, f_on, infer_div_on=true)          # Used as initial parameters in optimisation
        log_likelihood_para_3(para) = -log_likelihood(mc_p, para[1], Nf_p, mc_s, Nf_s, para[2], para[3], f_on)
        res = Optim.optimize(log_likelihood_para_3, [mu_off, mu_on, rel_div_on])                  # 3 inference parameters
        if Optim.converged(res) == true
            return [Optim.minimizer(res); 6 + 2*Optim.minimum(res)]                               # Mutation rate response-off/-on cells, relative division rate response-on cells, AIC
        else
            return [0., 0., -1., Inf]                                                                               
        end
    else
        mu_on = estimate_mu_het(mc_s, mu_off, Nf_s, f_on)[1]                                      # Used as initial parameter in optimisation
        log_likelihood_para_2(para) = -log_likelihood(mc_p, para[1], Nf_p, mc_s, Nf_s, para[2], rel_div_on, f_on)
        res = Optim.optimize(log_likelihood_para_2, [mu_off, mu_on])                              # 2 inference parameters
        if Optim.converged(res) == true
            return [Optim.minimizer(res); 4 + 2*Optim.minimum(res)]                               # Mutation rate response-off/-on cells, AIC
        else
            return [0., 0., Inf]                                                                                    
        end
    end
end


# The following inference functions are used to determine the initial parameters for the joint inference
function estimate_mu_hom(mc::Vector{Int}, Nf; fit_m=1.)                            # Estimating the mutation rate for a homogeneous population with optional differential fitness of mutants
    if fit_m == "infer"                                                            # Differential fitness of mutants is inferred 
        log_likelihood_para_2(para) = -log_likelihood(mc, para[1], para[2], Nf)
        res = Optim.optimize(log_likelihood_para_2, [initial_mu(mc, Nf, 100), 1.]) # 2 inference parameters
        if Optim.converged(res) == true
            return [Optim.minimizer(res); 4 + 2*Optim.minimum(res)]                # Mutation rate, differential fitness of mutants, AIC
        else
            return [0., -1., Inf]
        end                                                   
    else                                                                           # Differential fitness of mutants is given (default set to 1)
        log_likelihood_para_1(para) = -log_likelihood(mc, para, fit_m, Nf)
        res = Optim.optimize(log_likelihood_para_1, 0., maximum(mc)/Nf)            # 1 inference parameter
        if Optim.converged(res) == true
            return [Optim.minimizer(res), 2 + 2*Optim.minimum(res)]                # Mutation rate, AIC
        else
            return [0., Inf]                                                       
        end                                                                              
    end
end
function estimate_mu_het(mc::Vector{Int}, mu_off_per_div, Nf; infer_div_on=false)                                  # Estimation of the mutation rate of response-on cells for given mutation rate of reswponse-off cells and unknown fraction of response-on subpopulation
    if infer_div_on == false                                                                                       # Relative division rate of response-on cells is not inferred and set to 0 instead
        log_likelihood_para_2(para) = -log_likelihood(mc, mu_off_per_div, Nf, para[1], 0., para[2])
        res = Optim.optimize(log_likelihood_para_2, [initial_mu(mc, mu_off_per_div, Nf, 0.05, 100), 0.05])         # 2 inference parameters; initial value for the fraction of response-on subpopulation is set to 5%
        if Optim.converged(res) == true
            return Optim.minimizer(res)                                                                            # Mutation rate response-on cells, fraciton of response-on
        else
            return [0., 0.]
        end
    else                                                                                                           # Relative division rate of response-on cells is inferred
        log_likelihood_para_3(para) = -log_likelihood(mc, mu_off_per_div, Nf, para[1], para[2], para[3])
        res = Optim.optimize(log_likelihood_para_3, [initial_mu(mc, mu_off_per_div, Nf, 0.05, 100), 0.1, 0.05])    # 3 inference parameters; initial value for relative division rate of response-on cells is set to 0.1, for the fraction of response-on subpopulation is set to 5%
        if Optim.converged(res) == true
            return Optim.minimizer(res)                                                                            # Mutation rate of response-on cells, relative division rate response-on cells, fraction of response-on subpopulation
        else
            return [0., 0., 0.]
        end
    end
end
function estimate_mu_het(mc::Vector{Int}, mu_off_per_gen, Nf, f_on; infer_div_on=false)            # Estimation of the mutation rate of response-on cells for given mutation rate of response-off cells and known fraction of response-on subpopulation
    if infer_div_on == false                                                                       # Relative division rate of response-on cells is not inferred and set to 0 instead
        log_likelihood_para_1(para) = -log_likelihood(mc, mu_off_per_gen, Nf, para, 0., f_on)
        res = Optim.optimize(log_likelihood_para_1, 0., maximum(mc)/(Nf*f_on))                     # 1 inference parameter
        if Optim.converged(res) == true
            return Optim.minimizer(res)                                                            # Mutation rate response-on cells
        else
            return [0.]
        end
    else                                                                                                  # Relative division rate of response-on cells is inferred 
        log_likelihood_para_2(para) = -log_likelihood(mc, mu_off_per_gen, Nf, para[1], para[2], f_on)
        res = Optim.optimize(log_likelihood_para_2, [initial_mu(mc, mu_off_per_gen, Nf, f_on, 100), 0.1]) # 2 inference parameters; initial value for relative division rate of response-on cells is set to 0.1
        if Optim.converged(res) == true
            return Optim.minimizer(res)                                                                   # Mutation rate, relative division rate of response-on cells
        else
            return [0., 0.]
        end
    end
end

# Log-likelihood functions
# Log-likelihood to observe a mutant count mc for a homogeneous population with differential fitness of mutants
function log_likelihood(mc::Vector{Int}, mu_per_gen, fit_m, Nf) 
    if mu_per_gen<=0. || fit_m<0. || Nf<=0.                           # Boundaries of the parameter regime
        return -Inf
    else
        p = P_mutant_count(maximum(mc), mu_per_gen, Nf, fit_m=fit_m)
        return sum(counts(mc, 0:maximum(mc)) .* log.(p))
    end
end
# Joint log-likelihood to observe mutant counts under permissive/stressful conditions with a joint differential fitness of mutants
function log_likelihood(mc_p::Vector{Int}, mu_p_per_gen, Nf_p, mc_s::Vector{Int}, Nf_s, mu_s_per_gen, fit_m) 
    if mu_p_per_gen<=0. || Nf_p<=0. || Nf_s<=0. || mu_s_per_gen<=0. || fit_m<0.                                   # Boundaries of the parameter regime
        return -Inf
    else
        p_p = P_mutant_count(maximum(mc_p), mu_p_per_gen, Nf_p, fit_m=fit_m)                                      # Mutant count distribution without stress
        p_s = P_mutant_count(maximum(mc_s), mu_s_per_gen, Nf_s, fit_m=fit_m) 
        return sum(counts(mc_p, 0:maximum(mc_p)) .* log.(p_p)) + sum(counts(mc_s, 0:maximum(mc_s)) .* log.(p_s))  # The two observations are independent and their probabilities can be multiplied
    end
end
# Log-likelihood to observe a mutant count mc for a heterogeneous population with relative division rate of response-on cells 
function log_likelihood(mc::Vector{Int}, mu_off_per_div, Nf, mu_on_per_div, rel_div_on, f_on) 
    if mu_off_per_div<=0. || Nf<=0. || mu_on_per_div<=0. || rel_div_on<0. || f_on<0. || f_on>1.          # Boundaries of the parameter regime
        return -Inf
    else
        p = P_mutant_count(maximum(mc), mu_off_per_div, Nf, mu_on_per_div, f_on, rel_div_on=rel_div_on)
        return sum(counts(mc, 0:maximum(mc)) .* log.(p))
    end
end
# Joint log-likelihood to observe the mutant counts (a) mc_p and (b) mc_s for 
# (a) homogeneous population without differential fitness of mutants
# (b) heterogeneous population with relative division rate of response-on cells 
function log_likelihood(mc_p::Vector{Int}, mu_off_per_div, Nf_p, mc_s::Vector{Int}, Nf_s, mu_on_per_div, rel_div_on, f_on)  
    if mu_off_per_div<=0. || Nf_p<=0. || mu_on_per_div<=0. || rel_div_on<0. || f_on<0. || Nf_s<=0. || f_on>1.               # Boundaries of the parameter regime
        return -Inf
    else
        p_p = P_mutant_count(maximum(mc_p), mu_off_per_div, Nf_p)                                                           # Mutant count distribution: permissive
        p_s = P_mutant_count(maximum(mc_s), mu_off_per_div, Nf_s, mu_on_per_div, f_on, rel_div_on=rel_div_on)               # Mutant count distribution: stress
        return sum(counts(mc_p, 0:maximum(mc_p)) .* log.(p_p)) + sum(counts(mc_s, 0:maximum(mc_s)) .* log.(p_s))            # The two observations are independent and their probabilities can be multiplied
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
function initial_mu(z, mc::Vector{Int}, mu_off_per_div, Nf, f_on)             # Estimate the mutation rate under the heterogeneous-response model for given z   
    if z == 0.
        return (log(empirical_pgf(z, mc)) - mu_off_per_div*Nf*(1-f_on)) / (Nf*f_on)
    else
        return (log(empirical_pgf(z, mc))/(z-1) + mu_off_per_div*Nf*(1-f_on) * log(1-z)/z) / (Nf*f_on)
    end
end
function initial_mu(mc::Vector{Int}, mu_off_per_div, Nf, f_on, z_values::Int) # Estimate the mutation rate under the heterogeneous-response model by averaging over a number of z values 
    mu_on = 0.
    for i = 0:z_values-1
        mu_on += initial_mu(i/z_values, mc, mu_off_per_div, Nf, f_on)
    end
    return maximum([mu_on/z_values, 0.])
end