using SpecialFunctions, Optim, StatsBase, Roots

# The following recursive formulas for the pdf of the mutant count distribution are based on 
# Keller, P., & Antal, T. (2015). Mutant number distribution in an exponentially growing population. Journal of Statistical Mechanics: Theory and Experiment, 2015(1), P01011. https://doi.org/10.1088/1742-5468/2015/01/P01011

# Muntant count distribution for standard model with optional differential mutant fitness
# m: Number of mutations during growth phase
# fit_m: Mutant fitness (relative to non-mutants); default = 1
# Returns probabilities to observe 0,...,K mutants
function pdf_mudi(K::Int, m, fit_m=1.) 
    p = zeros(Float64, K+1)
    if fit_m == 0.
        for k = 0:K
            p[k+1] = m^k * exp(-m) / factorial(big(k))          # Fitness of mutants = 0 -> Poisson distribution
        end
    else
        q = Q(K, m, fit_m)
        p[1] = exp(q[1])
        for k = 1:K
            S = 0.
            for i = 0:k-1
                S += (k-i) * q[k-i+1] * p[i+1]
            end
            if isnan(S)                                         # For very small fit_m; approximated by fit_m=0 
                for k = 0:K
                    p[k+1] = m^k * exp(-m) / factorial(big(k))          
                end
                break
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

# Mutant count distribution for heterogeneous-response model (two subpopulations with stress response switched off/on, called off/on-cells). 
# m: Number of mutations arising in off-cells
# mu_het: Mutation-rate heterogeneity = mu_on/mu_off * f_on/(1-f_on)
# f_on: Fraction of on-cells at the end of the growth phase
# rel_div_on: Relative division rate of on-cells compared to off-cells
# Mutant count distribution not explicitely dependent on f_on if rel_div_on=0 (with the exception of the case f_on=0)
# Returns probabilities to observe 0,...,K mutants
function pdf_mudi(K::Int, m, mu_het, f_on, rel_div_on)
    if f_on == 0.
        return pdf_mudi(K, m)
    else
        if rel_div_on == 0.
            p_off = pdf_mudi(K, m)              # Contribution off-cells
            p_on = pdf_mudi(K, m * mu_het, 0.)  # Contribution on-cells
        else
            p_off = pdf_mudi(K, m * (1-f_on)/(1-f_on*(1-rel_div_on)))                                               # Contribution off-cells       
            p_on = pdf_mudi(K, m * mu_het * (1-f_on)/(1-f_on*(1-rel_div_on)), rel_div_on/(1-f_on*(1-rel_div_on)))   # Contribution on-cells
        end
        p = zeros(Float64, length(p_off))
        for k = 0:length(p_off)-1
            pk = 0.
            for j = 0:k
                pk += p_off[j+1] * p_on[k-j+1] # Total mutant count = sum of contributions from on-cells and the rest of the population
            end
            p[k+1] = pk
        end
        return p
    end
end
pdf_mudi(K::Int, m, mu_het, zero_div::Bool) = pdf_mudi(K, m, mu_het, 0.5, 0.)

# Mutation rate estimation algorithms 
# Return maximum likelihood estimates, confidence intervals, and AIC/BIC value (if the optimisation fails, Inf is returned)

# Mutation rate estimation from single fluctuation assay using the standard model with optional differential mutant fitness
# Input
# mc: Mutant counts
# Nf: Final population size
# Optional
# fit_m: Mutant fitness
#        By default fit_m=1 but can be set to different value if known from separate experiment
#        Alternatively, mutant fitness is inferred if fit_m=false
function estimu_hom(mc::Vector{Int}, Nf, fit_m::Float64=1.; conf=false)
	est_res = zeros(Float64, 9)
    est_res[4:6] = [fit_m, fit_m, fit_m]                                                                  
	log_likelihood_para_1(para) = -log_likelihood(mc, para, fit_m)                   # 1 inference parameter: Number of mutations
	res = Optim.optimize(log_likelihood_para_1, 0., maximum(mc))
	if Optim.converged(res) == true
        est_res[1] = Optim.minimizer(res)/Nf
        est_res[7] = Optim.minimum(res)
        est_res[8] = 2*Optim.minimum(res) + 2
        est_res[9] = 2*Optim.minimum(res) + log(length(mc))
        if conf
            b = CI(mc, Optim.minimizer(res), fit_m, Optim.minimum(res))    
		    est_res[2:3] = b./Nf
        end
	else                                                      
		est_res[7:9] = [Inf,Inf,Inf]
	end 
	return est_res
end 
function estimu_hom(mc::Vector{Int}, Nf, fit_m::Bool; conf=false)                                   # Mutant fitness not given -> inferred 
	est_res = zeros(Float64, 9)                                                         
	log_likelihood_para_2(para) = -log_likelihood(mc, para[1], para[2])                 # 2 inference parameters: Number of mutations, mutant fitness
	res = Optim.optimize(log_likelihood_para_2, [initial_m(mc, 1000), 1.]) 
	if Optim.converged(res) == true
		p = Optim.minimizer(res)
        est_res[1] = p[1]/Nf
        est_res[4] = p[2]
        est_res[7] = Optim.minimum(res)
        est_res[8] = 2*Optim.minimum(res) + 4
        est_res[9] = 2*Optim.minimum(res) + 2*log(length(mc))                   
        if conf
            b = CI(mc, p[1], p[2], true, Optim.minimum(res))
            est_res[2:3] = b[1,:]./Nf
            est_res[5:6] = b[2,:]
        end
	else
		est_res[7:9] = [Inf,Inf,Inf]
	end
	return est_res                                                         
end 

# Mutation rate estimation from pair of fluctuation assays under permissive/stressful cond. using the homogeneous-response model
# mc_p: Mutant counts under permissive cond.
# Nf_p: Average final population size under permissive cond.
# mc_s: Mutant counts under stressful cond.
# Nf_s: Average final population size under stressful cond.
# Optional
# fit_m: Mutant fitness
#        By default fit_m=1 but can be set to different value(s) if known from separate experiment(s)
#        Alternatively, can be set as inference parameter(s) via fit_m=false
#        If only one value is given, mutant fitness is constrained to be equal under permissive/stressful cond.
#        To not constrain mutant fitness, values have to be given as a tuple (mutant fitness permissive cond., mutant fitness stressful cond.)
function estimu_hom(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, fit_m::Float64=1.; conf=false)    # Mutant fitness set to input
    est_res = zeros(Float64, 21)
	est_res_p = estimu_hom(mc_p, Nf_p, fit_m, conf=conf)     # Estimation for permissive cond.
	est_res_s = estimu_hom(mc_s, Nf_s, fit_m, conf=conf)     # Estimation for stressful cond.
    est_res[1:12] = [est_res_p[1:6]; est_res_s[1:6]]
	if est_res_p[end] != Inf && est_res_s[end] != Inf
        est_res[13:15] = est_res_s[1:3] ./ [est_res_p[1], est_res_p[3], est_res_p[2]]
        est_res[16:18] = est_res_s[4:6] ./ [est_res_p[4], est_res_p[6], est_res_p[5]]
	    est_res[19:20] = est_res_p[7:8] .+ est_res_s[7:8]
        est_res[21] = 2*(est_res_p[7]+est_res_s[7]) + 2*log(length(mc_p)+length(mc_s))
    else
        est_res[19:21] = [Inf,Inf,Inf]
    end                                      
	return est_res
end
function estimu_hom(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, fit_m::Bool; conf=false)          # Mutant fitness not given -> inferred (as constrained parameter)
	est_res = zeros(Float64, 21)
	m_p, rho = estimate_init_hom(mc_p, false)                                                   
	m_s = estimate_init_hom(mc_s)                                                                    
	log_likelihood_para_3(para) = -log_likelihood(mc_p, para[1], mc_s, para[2], para[3])    # 3 inference parameters: Number of mutations under permissive/stressful cond., mutant fitness         
	res = Optim.optimize(log_likelihood_para_3, [m_p, m_s, rho])                                      
	if Optim.converged(res) == true
		p = Optim.minimizer(res)
        est_res[1] = p[1]/Nf_p
        est_res[4] = p[3]
        est_res[7] = p[2]/Nf_s
        est_res[10] = p[3]
        est_res[13] = p[2]/p[1]*Nf_p/Nf_s
        est_res[16:18] = [1.,1.,1.]
        est_res[19] = Optim.minimum(res)
        est_res[20] = 2*Optim.minimum(res) + 6
        est_res[21] = 2*Optim.minimum(res) + 3*log(length(mc_p)+length(mc_s))    
        if conf
            b = CI(mc_p, p[1], mc_s, p[2], p[3], true, Optim.minimum(res))  
            est_res[2:3] = b[1,:]./Nf_p
            est_res[5:6] = b[3,:]
            est_res[8:9] = b[2,:]./Nf_s  
            est_res[11:12] = b[3,:]
            est_res[14:15] = b[4,:].*(Nf_p/Nf_s)
        end   
	else
		est_res[19:21] = [Inf,Inf,Inf]
	end
	return est_res
end
function estimu_hom(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, fit_m::Tuple{<:Number,<:Number}; conf=false) # Mutant fitness not constrained 
    if typeof(fit_m[1]) == Float64 && fit_m[1] == fit_m[2]
		return estimu_hom(mc_p, Nf_p, mc_s, Nf_s, fit_m[1], conf=conf)
	else
        est_res = Vector{Float64}(undef, 21)
		fit_m_p = fit_m[1]
        fit_m_s = fit_m[2]
		est_res_p = estimu_hom(mc_p, Nf_p, fit_m_p, conf=conf)             # Estimation permissive cond.
		est_res_s = estimu_hom(mc_s, Nf_s, fit_m_s, conf=conf)             # Estimation stressful cond.
        est_res[1:12] = [est_res_p[1:6]; est_res_s[1:6]]
        if est_res_p[end] != Inf && est_res_s[end] != Inf       # Estimation under both conditions successful
            est_res[13:15] = est_res_s[1:3] ./ [est_res_p[1], est_res_p[3], est_res_p[2]]
            est_res[16:18] = est_res_s[4:6] ./ [est_res_p[4], est_res_p[6], est_res_p[5]]
            est_res[19:20] = est_res_p[7:8] .+ est_res_s[7:8]
            est_res[21] = 2*(est_res_p[7]+est_res_s[7]) + (2+sum([typeof(fit_m_p),typeof(fit_m_s)].==Bool))*log(length(mc_p)+length(mc_s))
        elseif est_res_p[end] != Inf && typeof(fit_m_s) == Bool # Estimation under stressful cond. failed -> set mutant fitness under stressful cond. to one
            est_res_s = estimu_hom(mc_s, Nf_s, conf=conf)
            if est_res_s[end] != Inf
            est_res[7:12] = est_res_s[1:6]
            est_res[13:15] = est_res_s[1:3] ./ [est_res_p[1], est_res_p[3], est_res_p[2]]
            est_res[16:18] = est_res_s[4:6] ./ [est_res_p[4], est_res_p[6], est_res_p[5]]
            est_res[19] = est_res_p[7] + est_res_s[7] 
            est_res[20] = est_res_p[8] + est_res_s[8] + 2
            est_res[21] = 2*(est_res_p[7]+est_res_s[7]) + (3+(typeof(fit_m_p)==Bool))*log(length(mc_p)+length(mc_s))
            end
        elseif est_res_s[end] != Inf && typeof(fit_m_p) == Bool # Estimation under permissive cond. failed -> set mutant fitness under permissive cond. to one
            est_res_p = estimu_hom(mc_p, Nf_p, conf=conf)
            if est_res_p[end] != Inf 
                est_res[1:6] = est_res_p[1:6]
                est_res[13:15] = est_res_s[1:3] ./ [est_res_p[1], est_res_p[3], est_res_p[2]]
                est_res[16:18] = est_res_s[4:6] ./ [est_res_p[4], est_res_p[6], est_res_p[5]]
                est_res[19:20] = est_res_p[7:8] .+ est_res_s[7:8]
                est_res[21] = 2*(est_res_p[7]+est_res_s[7]) + (3+(typeof(fit_m_s)==Bool))*log(length(mc_p)+length(mc_s))
            end
        else
            est_res[19:21] = [Inf,Inf,Inf]
        end                      
		return est_res
	end
end

# Mutation rate estimation from pair of fluctuation assays under permissive/stressful cond. using the heterogeneous-response model 
# mc_p: Mutant counts under permissive cond.
# Nf_p: Average final population size under permissive cond.
# mc_s: Mutant counts under stressful cond.
# Nf_s: Average final population size under stressful cond.
# Optional
# f_on: Fraction of on-cells
#       By default inferred, can be set to different value if known from separate experiment
# rel_div_on: Relative division rate of on-cells compared to off-cells
#             By default set to rel_div_on=0., inferred if rel_div_on=false
#             For rel_div_on=0, the fraction of on-cells cannot be inferred
function estimu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, f_on::Float64, rel_div_on::Float64=0.; conf=false) # Fraction and relative division rate of on-cells given
    est_res = zeros(Float64, 24)	
	m = estimate_init_hom(mc_p)*Nf_s/Nf_p                       # Initial value for optimisation   
    mu_het, div_init = estimate_init_het(mc_s, m, f_on, rel_div_on)                       
    log_likelihood_para_2(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], f_on, rel_div_on)    # 2 inference parameters: Number of mutations in off-cells, mutation-rate heterogeneity
    res = Optim.optimize(log_likelihood_para_2, [m, mu_het])                                                 
    if Optim.converged(res) == true
        p = Optim.minimizer(res)
        est_res[1] = p[1]/Nf_s
        est_res[4] = p[2]
        est_res[7:9] = [f_on,f_on,f_on]
        est_res[10:12] = [rel_div_on,rel_div_on,rel_div_on]
        est_res[13] = p[2]*p[1]*(1-f_on)/(f_on*Nf_s)
        est_res[16] = p[2]*(1-f_on)/f_on
        est_res[19] = (1-f_on)*(1+p[2])
        est_res[22] = Optim.minimum(res) 
        est_res[23] = 2*Optim.minimum(res) + 4
        est_res[24] = 2*Optim.minimum(res) + 2*log(length(mc_p)+length(mc_s))
        if conf
            b = CI(mc_p, mc_s, Nf_p/Nf_s, p[1], p[2], f_on, rel_div_on, Optim.minimum(res))
            est_res[2:3] = b[1,:]./Nf_s
            est_res[5:6] = b[2,:]
            est_res[14:15] = b[3,:]./Nf_s
            est_res[17:18] = b[2,:].*((1-f_on)/f_on)
            est_res[20:21] = (1-f_on).*(1 .+b[2,:])
        end                                                        
    else
        est_res[22:24] = [Inf,Inf,Inf]
    end   
    return est_res 
end
function estimu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, f_on::Float64, rel_div_on::Bool; conf=false)    # Relative division rate on-cells not given -> inferred
    est_res = zeros(Float64, 24)	
	m = estimate_init_hom(mc_p)*Nf_s/Nf_p  
    mu_het, rel_div_on = estimate_init_het(mc_s, m, f_on)
    log_likelihood_para_3(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], f_on, para[3])     
    res = Optim.optimize(log_likelihood_para_3, [m, mu_het, rel_div_on])                                    # 3 inference parameters: Number of mutations in off-cells, mutation-rate heterogeneity, relative division rate on-cells                         
    if Optim.converged(res) == true
        p = Optim.minimizer(res)
        est_res[1] = p[1]/Nf_s
        est_res[4] = p[2]
        est_res[7:9] = [f_on,f_on,f_on]
        est_res[10] = p[3]
        est_res[13] = p[2]*p[1]*(1-f_on)/(f_on*Nf_s)
        est_res[16] = p[2]*(1-f_on)/f_on
        est_res[19] = (1-f_on)*(1+p[2])
        est_res[22] = Optim.minimum(res) 
        est_res[23] = 2*Optim.minimum(res) + 6
        est_res[24] = 2*Optim.minimum(res) + 3*log(length(mc_p)+length(mc_s))
        if conf
            b = CI(mc_p, mc_s, Nf_p/Nf_s, p[1], p[2], f_on, p[3], true, Optim.minimum(res))
            est_res[2:3] = b[1,:]./Nf_s
            est_res[5:6] = b[2,:]
            est_res[11:12] = b[3,:]
            est_res[14:15] = b[4,:]./Nf_s
            est_res[17:18] = b[2,:].*((1-f_on)/f_on)
            est_res[20:21] = (1-f_on).*(1 .+b[2,:])
        end
    else
        est_res[22:24] = [Inf,Inf,Inf]
    end
    return est_res
end
function estimu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, f_on::Bool, rel_div_on::Float64=0.; conf=false)   # Fraction of on-cells not given
    m = estimate_init_hom(mc_p)*Nf_s/Nf_p 
    if rel_div_on == 0.                                                                                     # Zero relative division rate on-cells -> Fraction of on-cells cannot be inferred   
        est_res = zeros(Float64, 9)
        log_likelihood_para_2(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2])     
        res = Optim.optimize(log_likelihood_para_2, [m, initial_mu_het(mc_s, m, 1000)])                     # 2 inference parameters: Number of mutations in off-cells, mutation-rate heterogeneity
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            est_res[1] = p[1]/Nf_s
            est_res[4] = p[2]
            est_res[7] = Optim.minimum(res)
            est_res[8] = 2*Optim.minimum(res) + 4
            est_res[9] = 2*Optim.minimum(res) + 2*log(length(mc_p)+length(mc_s))
            if conf
                b = CI(mc_p, mc_s, Nf_p/Nf_s, p[1], p[2], Optim.minimum(res))
                est_res[2:3] = b[1,:]./Nf_s
                est_res[5:6] = b[2,:]
            end                                                            
        else
            est_res[7:9] = [Inf,Inf,Inf]
        end
    else                                                                                                    # Non-zero relative division rate on-cells -> Fraction of on-cells inferred
        est_res = zeros(Float64, 24)
        mu_het, div_init = estimate_init_het(mc_s, m, 0., rel_div_on)                                   
        mu_inc = estimate_init_hom(mc_s)/m                                                                       
        f_upper = 1 - mu_inc/(mu_het+1)                                                          
        if f_upper <= 0.
            f_upper = 1/mu_inc
        end
        f_lower = - log(1-f_upper) / log(Nf_s)                              # Initial value for optimisation                                                                                                                                                     
        mu_het, div_init = estimate_init_het(mc_s, m, f_lower, div_init)                                        
        f_on = f_lower/(1 - rel_div_on)                                                                         
        f_on = minimum([maximum([f_lower, f_on]), f_upper])  
        log_likelihood_para_3(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], para[3], rel_div_on)     
        res = Optim.optimize(log_likelihood_para_3, [m, mu_het, f_on])                                              # 3 inference parameters: Number of mutations in off-cells, mutation-rate heterogeneity, fraction of on-cells                         
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            est_res[1] = p[1]/Nf_s
            est_res[4] = p[2]
            est_res[7] = p[3]
            est_res[10:12] = [rel_div_on,rel_div_on,rel_div_on]
            est_res[13] = p[2]*p[1]*(1-p[3])/(p[3]*Nf_s)
            est_res[16] = p[2]*(1-p[3])/p[3]
            est_res[19] = (1-p[3])*(1+p[2])
            est_res[22] = Optim.minimum(res) 
            est_res[23] = 2*Optim.minimum(res) + 6
            est_res[24] = 2*Optim.minimum(res) + 3*log(length(mc_p)+length(mc_s))
            if conf
                b = CI(mc_p, mc_s, Nf_p/Nf_s, p[1], p[2], p[3], rel_div_on, false, Optim.minimum(res))
                est_res[2:3] = b[1,:]./Nf_s
                est_res[5:6] = b[2,:]
                est_res[8:9] = b[3,:]
                est_res[14:15] = b[4,:]./Nf_s
                est_res[17:18] = b[5,:]
                est_res[20:21] = b[6,:]
            end                                                           
        else
            est_res[22:24] = [Inf,Inf,Inf]
        end   
    end
    return est_res 
end
function estimu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, f_on::Bool, rel_div_on::Bool; conf=false)    # Relative division rate on-cells and fraction of on-cells not given -> inferred
    est_res = zeros(Float64, 24)
    m = estimate_init_hom(mc_p)*Nf_s/Nf_p    
    mu_het, rel_div_on = estimate_init_het(mc_s, m, 0.)                                   
    mu_inc = estimate_init_hom(mc_s)/m                                                                       
    f_upper = 1 - mu_inc/(mu_het+1)                                                     
    if f_upper <= 0.
        f_upper = 1/mu_inc
    end
    f_lower = - log(1-f_upper) / log(Nf_s)                                                                                                                                                           
    mu_het, rel_div_on = estimate_init_het(mc_s, m, f_lower, rel_div_on)                                        
    f_on = f_lower/(1 - rel_div_on)                                                                         
    f_on = minimum([maximum([f_lower, f_on]), f_upper])  
    log_likelihood_para_4(para) = -log_likelihood(mc_p, mc_s, Nf_p/Nf_s, para[1], para[2], para[3], para[4])     
    res = Optim.optimize(log_likelihood_para_4, [m, mu_het, f_on, rel_div_on])                                    # 4 inference parameters: Number of mutations in off-cells, mutation-rate heterogeneity, relative division rate on-cells, fraction of on-cells                         
    if Optim.converged(res) == true
        p = Optim.minimizer(res)
        est_res[1] = p[1]/Nf_s
        est_res[4] = p[2]
        est_res[7] = p[3]
        est_res[10] = p[4]
        est_res[13] = p[2]*p[1]*(1-p[3])/(p[3]*Nf_s)
        est_res[16] = p[2]*(1-p[3])/p[3]
        est_res[19] = (1-p[3])*(1+p[2])
        est_res[22] = Optim.minimum(res) 
        est_res[23] = 2*Optim.minimum(res) + 8
        est_res[24] = 2*Optim.minimum(res) + 4*log(length(mc_p)+length(mc_s))
        if conf
            b = CI(mc_p, mc_s, Nf_p/Nf_s, p[1], p[2], p[3], p[4], true, true, Optim.minimum(res))
            est_res[2:3] = b[1,:]./Nf_s
            est_res[5:6] = b[2,:]
            est_res[8:9] = b[3,:]
            est_res[11:12] = b[4,:]
            est_res[14:15] = b[5,:]./Nf_s
            est_res[17:18] = b[6,:]
            est_res[20:21] = b[7,:]
        end
    else
        est_res[22:24] = [Inf,Inf,Inf]
    end
    return est_res
end

# Inference functions to calculate initial parameters for the joint ML estimation
# Estimating the number of mutations for a homogeneous population with mutant fitness fit_m (inferred if = false)
function estimate_init_hom(mc::Vector{Int}, fit_m::Float64=1.)              # Mutant fitness given (default = 1)
    log_likelihood_para_1(para) = -log_likelihood(mc, para, fit_m)          # 1 inference parameter: Number of mutations
    res = Optim.optimize(log_likelihood_para_1, 0., maximum(mc))
    if Optim.converged(res) == true
        return Optim.minimizer(res)                                        
    else
        return 1.                                                       
    end                                                                              
end
function estimate_init_hom(mc::Vector{Int}, fit_m::Bool)                    # Mutant fitness not given -> inferred                                                                          
    log_likelihood_para_2(para) = -log_likelihood(mc, para[1], para[2])     # 2 inference parameters: Number of mutations, mutant fitness
    res = Optim.optimize(log_likelihood_para_2, [initial_m(mc, 1000), 1.]) 
    if Optim.converged(res) == true
        return Optim.minimizer(res)                                        
    else
        return [estimate_init_hom(mc), 1.]
    end                                                                                                                                
end
# Estimating the mutation-rate heterogeneity for given number of mutations (under stressful cond.) and known fraction of on-cells
# Optional input parameter: Initial value for the relative division rate of on-cells
function estimate_init_het(mc::Vector{Int}, m, f_on, rel_div_on=0.)                                       
    log_likelihood_para_2(para) = -log_likelihood(mc, m, para[1], f_on, para[2])            # 2 inference parameters: Mutation-rate heterogeneity, relative division rate of on-cells                
    res = Optim.optimize(log_likelihood_para_2, [initial_mu_het(mc, m, 1000), rel_div_on])          
    if Optim.converged(res) == true
        return Optim.minimizer(res)                                                                                                                    
    else
        return [maximum([initial_mu_het(mc, m, 1000), f_on]), rel_div_on]                   # If inference fails, input parameter for relative division rate of on-cells is returned
    end
end

# Log-likelihood functions
# Log-likelihood to observe mutant counts mc for standard model with mutant fitness
function log_likelihood(mc::Vector{Int}, m, fit_m) 
    if m<=0. || fit_m<0.                            # Boundaries of the parameter regime
        return -Inf
    else
        p = pdf_mudi(maximum(mc), m, fit_m)
        ll = sum(counts(mc, 0:maximum(mc)) .* log.(p))
        if !isnan(ll)
            return ll
        else 
            return -Inf
        end
    end
end
# Joint log-likelihood to observe mutant counts for standard model with mutant fitness constrained to be equal under permissive/stressful conditions
# m_p: Number of mutations during growth phase under permissive cond.
# m_s: Number of mutations during growth phase under stressful cond.
function log_likelihood(mc_p::Vector{Int}, m_p, mc_s::Vector{Int}, m_s, fit_m) 
    if m_p<=0. || m_s<=0. || fit_m<0.               # Boundaries of the parameter regime
        return -Inf
    else
        p_p = pdf_mudi(maximum(mc_p), m_p, fit_m)   # Under permissive cond.                                               
        p_s = pdf_mudi(maximum(mc_s), m_s, fit_m)   # Under stressful cond.
        ll = sum(counts(mc_p, 0:maximum(mc_p)) .* log.(p_p)) + sum(counts(mc_s, 0:maximum(mc_s)) .* log.(p_s)) 
        if !isnan(ll)
            return ll
        else 
            return -Inf
        end
    end
end
# Log-likelihood to observe mutant counts mc for a heterogeneous population
function log_likelihood(mc::Vector{Int}, m, mu_het, f_on, rel_div_on) 
    if m<=0. || mu_het<0. || rel_div_on<0. || f_on<0. || f_on>=1.           # Boundaries of the parameter regime
        return -Inf
    else
        p = pdf_mudi(maximum(mc), m, mu_het, f_on, rel_div_on)
        ll = sum(counts(mc, 0:maximum(mc)) .* log.(p))
        if !isnan(ll)
            return ll
        else 
            return -Inf
        end
    end
end
# Joint log-likelihood to observe mutant counts for heterogeneous-response model 
# N_ratio: Ratio of final population size permissive cond. over final population size stressful cond.
function log_likelihood(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio, m, mu_het, f_on, rel_div_on)  
    if m<=0. || mu_het<0. || rel_div_on<0. || f_on<0. || f_on>=1.           # Boundaries of the parameter regime
        return -Inf
    else
        p_p = pdf_mudi(maximum(mc_p), m*N_ratio)                            # Under permissive cond.
        p_s = pdf_mudi(maximum(mc_s), m, mu_het, f_on, rel_div_on)          # Under stressful cond.
        ll = sum(counts(mc_p, 0:maximum(mc_p)) .* log.(p_p)) + sum(counts(mc_s, 0:maximum(mc_s)) .* log.(p_s))    
        if !isnan(ll)
            return ll
        else
            return -Inf
        end
    end
end
# Joint log-likelihood to observe mutant counts for heterogeneous-response model with zero relative division rate of on-cells
function log_likelihood(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio, m, mu_het)  
    if m<=0. || mu_het<0.                                       # Boundaries of the parameter regime
        return -Inf
    else
        p_p = pdf_mudi(maximum(mc_p), m*N_ratio)                # Under permissive cond.
        p_s = pdf_mudi(maximum(mc_s), m, mu_het, true)          # Under stressful cond.
        ll = sum(counts(mc_p, 0:maximum(mc_p)) .* log.(p_p)) + sum(counts(mc_s, 0:maximum(mc_s)) .* log.(p_s))    
        if !isnan(ll)
            return ll
        else
            return -Inf
        end
    end
end

# Profile log-likelihood to calculate confidence intervals
# The 95% quantile of the Chi-squared distribution used to calculate CIs
# Depend on the observed mutant counts, on which parameters are estimated, the ML estimates and the maximum likelihood itself 
chisq_95 = 3.84145882069412447634704221854917705059051513671875
function CI(mc::Vector{Int}, m, fit_m, ML)
    function log_likelihood_ratio_1(para) 
        if para == m
            return -0.5*chisq_95
        else
            return -log_likelihood(mc, para, fit_m) - ML - 0.5*chisq_95
        end
    end
    l_1 = find_zero(log_likelihood_ratio_1, (0., m))
    u_1 = m
    try
        u_1 = find_zero(log_likelihood_ratio_1, (m, maximum(mc)))
    catch err
        u_1 = find_zero(log_likelihood_ratio_1, (m, 10*maximum(mc)))
    end
    return [l_1 u_1]
end
function CI(mc::Vector{Int}, m, fit_m, infer_fit_m::Bool, ML)
    function log_likelihood_ratio_1(para)
        if para == m
            return -0.5*chisq_95
        else
            log_likelihood_1(P) = -log_likelihood(mc, para, P[1])
            res = Optim.optimize(log_likelihood_1, [fit_m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood(mc, para, P1) - ML - 0.5*chisq_95
        end
    end
    l_1 = find_zero(log_likelihood_ratio_1, (0., m))
    u_1 = m
    try
        u_1 = find_zero(log_likelihood_ratio_1, (m, maximum(mc)))
    catch err
        u_1 = find_zero(log_likelihood_ratio_1, (m, 10*maximum(mc)))
    end
    function log_likelihood_ratio_2(para) 
        if para == fit_m
            return -0.5*chisq_95
        else
            log_likelihood_2(P) = -log_likelihood(mc, P[1], para)
            res = Optim.optimize(log_likelihood_2, [m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood(mc, P1, para) - ML - 0.5*chisq_95
        end
    end
    l_2 = 0.
    try
        l_2 = find_zero(log_likelihood_ratio_2, (0., fit_m))
    catch err
    end
    u_2 = find_zero(log_likelihood_ratio_2, (fit_m, 100.))
    return [l_1 u_1; l_2 u_2]
end
function CI(mc_p::Vector{Int}, m_p, mc_s::Vector{Int}, m_s, fit_m, ML)
    if length(fit_m) == 2
        fit_m_p, fit_m_s = fit_m
    else
        fit_m_p, fit_m_s = fit_m, fit_m
    end
    function log_likelihood_ratio_1(para) 
        if para == m_p
            return -0.5*chisq_95
        else
            return -log_likelihood(mc_p, para, fit_m_p) -log_likelihood(mc_s, m_s, fit_m_s) - ML - 0.5*chisq_95
        end
    end
    l_1 = find_zero(log_likelihood_ratio_1, (0., m_p))
    u_1 = m_p
    try
        u_1 = find_zero(log_likelihood_ratio_1, (m_p, maximum(mc_p)))
    catch err
        u_1 = find_zero(log_likelihood_ratio_1, (m_p, 10*maximum(mc_p)))
    end
    function log_likelihood_ratio_2(para)
        if para == m_s
            return -0.5*chisq_95 
        else
            return -log_likelihood(mc_p, m_p, fit_m_p) -log_likelihood(mc_s, para, fit_m_s) - ML - 0.5*chisq_95
        end
    end
    l_2 = find_zero(log_likelihood_ratio_2, (0., m_s))
    u_2 = m_s
    try
        u_2 = find_zero(log_likelihood_ratio_2, (m_s, maximum(mc_s)))
    catch err
        u_2 = find_zero(log_likelihood_ratio_2, (m_s, 10*maximum(mc_s)))
    end
    return [l_1 u_1; l_2 u_2]
end
function CI(mc_p::Vector{Int}, m_p, mc_s::Vector{Int}, m_s, fit_m, infer_fit_m::Bool, ML)
    del_mu = Vector{Float64}(undef, 6)
    function log_likelihood_ratio_1(para)
        if para == m_p
            return -0.5*chisq_95
        else
            log_likelihood_1(P) = -log_likelihood(mc_p, para, mc_s, P[1], P[2])
            res = Optim.optimize(log_likelihood_1, [m_s, fit_m])
            P1, P2 = Optim.minimizer(res)
            return -log_likelihood(mc_p, para, mc_s, P1, P2) - ML - 0.5*chisq_95
        end
    end
    l_1 = find_zero(log_likelihood_ratio_1, (0., m_p))
    l_1_del_mu(P) = -log_likelihood(mc_p, l_1, mc_s, P[1], P[2])
    res = Optim.optimize(l_1_del_mu, [m_s, fit_m])
    del_mu[1] = Optim.minimizer(res)[1]/l_1
    u_1 = m_p
    try
        u_1 = find_zero(log_likelihood_ratio_1, (m_p, maximum(mc_p)))
    catch err
        u_1 = find_zero(log_likelihood_ratio_1, (m_p, 10*maximum(mc_p)))
    end
    u_1_del_mu(P) = -log_likelihood(mc_p, u_1, mc_s, P[1], P[2])
    res = Optim.optimize(l_1_del_mu, [m_s, fit_m])
    del_mu[2] = Optim.minimizer(res)[1]/u_1
    function log_likelihood_ratio_2(para)
        if para == m_s
            return -0.5*chisq_95
        else
            log_likelihood_2(P) = -log_likelihood(mc_p, P[1], mc_s, para, P[2])
            res = Optim.optimize(log_likelihood_2, [m_p, fit_m])
            P1, P2 = Optim.minimizer(res)
            return -log_likelihood(mc_p, P1, mc_s, para, P2) - ML - 0.5*chisq_95
        end
    end
    l_2 = find_zero(log_likelihood_ratio_2, (0., m_s))
    l_2_del_mu(P) = -log_likelihood(mc_p, P[1], mc_s, l_2, P[2])
    res = Optim.optimize(l_2_del_mu, [m_p, fit_m])
    del_mu[3] = l_2/Optim.minimizer(res)[1]
    u_2 = m_s
    try
        u_2 = find_zero(log_likelihood_ratio_2, (m_s, maximum(mc_s)))
    catch err
        u_2 = find_zero(log_likelihood_ratio_2, (m_s, 10*maximum(mc_s)))
    end
    u_2_del_mu(P) = -log_likelihood(mc_p, P[1], mc_s, u_2, P[2])
    res = Optim.optimize(u_2_del_mu, [m_p, fit_m])
    del_mu[4] = u_2/Optim.minimizer(res)[1]
    function log_likelihood_ratio_3(para) 
        if para == fit_m
            return -0.5*chisq_95
        else
            log_likelihood_3(P) = -log_likelihood(mc_p, P[1], mc_s, P[2], para)
            res = Optim.optimize(log_likelihood_3, [m_p, m_s])
            P1, P2 = Optim.minimizer(res)
            return -log_likelihood(mc_p, P1, mc_s, P2, para) - ML - 0.5*chisq_95
        end
    end
    l_3 = 0.
    try
        l_3 = find_zero(log_likelihood_ratio_3, (0., fit_m))
    catch err
    end
    l_3_del_mu(P) = -log_likelihood(mc_p, P[1], mc_s, P[2], l_3)
    res = Optim.optimize(l_3_del_mu, [m_p, m_s])
    del_mu[5] = Optim.minimizer(res)[2]/Optim.minimizer(res)[1]
    u_3 = find_zero(log_likelihood_ratio_3, (fit_m, Inf))
    u_3_del_mu(P) = -log_likelihood(mc_p, P[1], mc_s, P[2], u_3)
    res = Optim.optimize(u_3_del_mu, [m_p, m_s])
    del_mu[6] = Optim.minimizer(res)[2]/Optim.minimizer(res)[1]
    return [l_1 u_1; l_2 u_2; l_3 u_3; minimum(del_mu) maximum(del_mu)]
end
function CI(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio, m, mu_het, f_on, rel_div_on, ML)
    m_on = Vector{Float64}(undef, 4)
    function log_likelihood_ratio_1(para)
        if para == m
            return -0.5*chisq_95
        else
            log_likelihood_1(P) = -log_likelihood(mc_p, mc_s, N_ratio, para, P[1], f_on, rel_div_on)
            res = Optim.optimize(log_likelihood_1, [mu_het])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood(mc_p, mc_s, N_ratio, para, P1, f_on, rel_div_on) - ML - 0.5*chisq_95
        end
    end
    l_1 = find_zero(log_likelihood_ratio_1, (0., m))
    l_1_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, l_1, P[1], f_on, rel_div_on)
    res = Optim.optimize(l_1_P, [mu_het])
    m_on[1] = Optim.minimizer(res)[1]*l_1*(1-f_on)/f_on
    u_1 = m
    try
        u_1 = find_zero(log_likelihood_ratio_1, (m, maximum([mc_p; mc_s])))
    catch err
        u_1 = find_zero(log_likelihood_ratio_1, (m, 10*maximum([mc_p; mc_s])))
    end
    u_1_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, u_1, P[1], f_on, rel_div_on)
    res = Optim.optimize(u_1_P, [mu_het])
    m_on[2] = Optim.minimizer(res)[1]*u_1*(1-f_on)/f_on
    function log_likelihood_ratio_2(para)
        if para == mu_het
            return -0.5*chisq_95
        else 
            log_likelihood_2(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], para, f_on, rel_div_on)
            res = Optim.optimize(log_likelihood_2, [m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood(mc_p, mc_s, N_ratio, P1, para, f_on, rel_div_on) - ML - 0.5*chisq_95
        end
    end
    l_2 = 0.
    try 
        l_2 = find_zero(log_likelihood_ratio_2, (0., mu_het))
    catch err
    end
    l_2_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, P, l_2, f_on, rel_div_on)
    res = Optim.optimize(l_2_P, 0., maximum([mc_p; mc_s]))
    m_on[3] = l_2*Optim.minimizer(res)[1]*(1-f_on)/f_on
    u_2 = find_zero(log_likelihood_ratio_2, (mu_het, Inf))
    u_2_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, P, u_2, f_on, rel_div_on)
    res = Optim.optimize(u_2_P, 0., maximum([mc_p; mc_s]))
    m_on[4] = u_2*Optim.minimizer(res)[1]*(1-f_on)/f_on
    return [l_1 u_1; l_2 u_2; minimum(m_on) maximum(m_on)]
end
function CI(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio, m, mu_het, f_on, rel_div_on, infer_r_f::Bool, ML)
    if infer_r_f == true
        m_on = Vector{Float64}(undef, 6)
        function log_likelihood_ratio_r1(para)
            if para == m
                return -0.5*chisq_95
            else
                log_likelihood_r1(P) = -log_likelihood(mc_p, mc_s, N_ratio, para, P[1], f_on, P[2])
                res = Optim.optimize(log_likelihood_r1, [mu_het, rel_div_on])
                P1, P2 = Optim.minimizer(res)
                return -log_likelihood(mc_p, mc_s, N_ratio, para, P1, f_on, P2) - ML - 0.5*chisq_95
            end
        end
        l_1 = find_zero(log_likelihood_ratio_r1, (0., m))
        l_1_rP(P) = -log_likelihood(mc_p, mc_s, N_ratio, l_1, P[1], f_on, P[2])
        res = Optim.optimize(l_1_rP, [mu_het, rel_div_on])
        m_on[1] = Optim.minimizer(res)[1]*l_1*(1-f_on)/f_on
        u_1 = m
        try
            u_1 = find_zero(log_likelihood_ratio_r1, (m, maximum([mc_p; mc_s])))
        catch err
            u_1 = find_zero(log_likelihood_ratio_r1, (m, 10*maximum([mc_p; mc_s])))
        end
        u_1_rP(P) = -log_likelihood(mc_p, mc_s, N_ratio, u_1, P[1], f_on, P[2])
        res = Optim.optimize(u_1_rP, [mu_het, rel_div_on])
        m_on[2] = Optim.minimizer(res)[1]*u_1*(1-f_on)/f_on
        function log_likelihood_ratio_r2(para) 
            if para == mu_het
                return -0.5*chisq_95
            else
                log_likelihood_r2(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], para, f_on, P[2])
                res = Optim.optimize(log_likelihood_r2, [m, rel_div_on])
                P1, P2 = Optim.minimizer(res)
                return -log_likelihood(mc_p, mc_s, N_ratio, P1, para, f_on, P2) - ML - 0.5*chisq_95
            end
        end
        l_2 = 0.
        try 
            l_2 = find_zero(log_likelihood_ratio_r2, (0., mu_het))
        catch err
        end
        l_2_rP(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], l_2, f_on, P[2])
        res = Optim.optimize(l_2_rP, [m, rel_div_on])
        m_on[3] = l_2*Optim.minimizer(res)[1]*(1-f_on)/f_on
        u_2 = find_zero(log_likelihood_ratio_r2, (mu_het, Inf))
        u_2_rP(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], u_2, f_on, P[2])
        res = Optim.optimize(u_2_rP, [m, rel_div_on])
        m_on[4] = u_2*Optim.minimizer(res)[1]*(1-f_on)/f_on
        function log_likelihood_ratio_r3(para)
            if para == rel_div_on
                return -0.5*chisq_95
            else 
                log_likelihood_r3(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], f_on, para)
                res = Optim.optimize(log_likelihood_r3, [m, mu_het])
                P1, P2 = Optim.minimizer(res)
                return -log_likelihood(mc_p, mc_s, N_ratio, P1, P2, f_on, para) - ML - 0.5*chisq_95
            end
        end
        l_3 = 0.
        try
            l_3 = find_zero(log_likelihood_ratio_r3, (0., rel_div_on))
        catch err
        end
        l_3_rP(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], f_on, l_3)
        res = Optim.optimize(l_3_rP, [m, mu_het])
        m_on[5] = Optim.minimizer(res)[2]*Optim.minimizer(res)[1]*(1-f_on)/f_on
        u_3 = find_zero(log_likelihood_ratio_r3, (rel_div_on, Inf))
        u_3_rP(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], f_on, u_3)
        res = Optim.optimize(u_3_rP, [m, mu_het])
        m_on[6] = Optim.minimizer(res)[2]*Optim.minimizer(res)[1]*(1-f_on)/f_on
        return [l_1 u_1; l_2 u_2; l_3 u_3; minimum(m_on) maximum(m_on)]
    else
        m_on = Vector{Float64}(undef, 6)
        mu_inc = Vector{Float64}(undef, 6)
        del_mu = Vector{Float64}(undef, 6)
        function log_likelihood_ratio_f1(para)
            if para == m
                return -0.5*chisq_95
            else
                log_likelihood_f1(P) = -log_likelihood(mc_p, mc_s, N_ratio, para, P[1], P[2], rel_div_on)
                res = Optim.optimize(log_likelihood_f1, [mu_het, f_on])
                P1, P2 = Optim.minimizer(res)
                return -log_likelihood(mc_p, mc_s, N_ratio, para, P1, P2, f_on) - ML - 0.5*chisq_95
            end
        end
        l_1 = 0.
        try
            l_1 = find_zero(log_likelihood_ratio_f1, (0., m))
        catch err
        end
        l_1_fP(P) = -log_likelihood(mc_p, mc_s, N_ratio, l_1, P[1], P[2], rel_div_on)
        res = Optim.optimize(l_1_fP, [mu_het, f_on])
        m_on[1] = Optim.minimizer(res)[1]*l_1*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
        mu_inc[1] = Optim.minimizer(res)[1]*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
        del_mu[1] = (1-Optim.minimizer(res)[2])*(1+Optim.minimizer(res)[1])
        u_1 = m
        try
            u_1 = find_zero(log_likelihood_ratio_f1, (m, maximum([mc_p; mc_s])))
        catch err
            u_1 = find_zero(log_likelihood_ratio_f1, (m, 10*maximum([mc_p; mc_s])))
        end
        u_1_fP(P) = -log_likelihood(mc_p, mc_s, N_ratio, u_1, P[1], P[2], rel_div_on)
        res = Optim.optimize(u_1_fP, [mu_het, f_on])
        m_on[2] = Optim.minimizer(res)[1]*u_1*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
        mu_inc[2] = Optim.minimizer(res)[1]*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
        del_mu[2] = (1-Optim.minimizer(res)[2])*(1+Optim.minimizer(res)[1])
        function log_likelihood_ratio_f2(para) 
            if para == mu_het
                return -0.5*chisq_95
            else
                log_likelihood_f2(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], para, P[2], rel_div_on)
                res = Optim.optimize(log_likelihood_f2, [m, f_on])
                P1, P2 = Optim.minimizer(res)
                return -log_likelihood(mc_p, mc_s, N_ratio, P1, para, P2, rel_div_on) - ML - 0.5*chisq_95
            end
        end
        l_2 = 0.
        try 
            l_2 = find_zero(log_likelihood_ratio_f2, (0., mu_het))
        catch err
        end
        l_2_fP(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], l_2, P[2], rel_div_on)
        res = Optim.optimize(l_2_fP, [m, f_on])
        m_on[3] = l_2*Optim.minimizer(res)[1]*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
        mu_inc[3] = l_2*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
        del_mu[3] = (1-Optim.minimizer(res)[2])*(1+l_2)
        u_2 = find_zero(log_likelihood_ratio_f2, (mu_het, Inf))
        u_2_fP(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], u_2, P[2], rel_div_on)
        res = Optim.optimize(u_2_fP, [m, f_on])
        m_on[4] = u_2*Optim.minimizer(res)[1]*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
        mu_inc[4] = u_2*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
        del_mu[4] = (1-Optim.minimizer(res)[2])*(1+u_2)
        function log_likelihood_ratio_f3(para)
            if para == rel_div_on
                return -0.5*chisq_95
            else 
                log_likelihood_f3(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], para, rel_div_on)
                res = Optim.optimize(log_likelihood_r3, [m, mu_het])
                P1, P2 = Optim.minimizer(res)
                return -log_likelihood(mc_p, mc_s, N_ratio, P1, P2, para, rel_div_on) - ML - 0.5*chisq_95
            end
        end
        l_3 = 0.
        try
            l_3 = find_zero(log_likelihood_ratio_f3, (0., f_on))
        catch err
        end
        l_3_fP(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], l_3, rel_div_on)
        res = Optim.optimize(l_3_fP, [m, mu_het])
        m_on[5] = Optim.minimizer(res)[1]*Optim.minimizer(res)[2]*(1-l_3)/l_3
        mu_inc[5] = Optim.minimizer(res)[2]*(1-l_3)/l_3
        del_mu[5] = (1-l_3)*(1+Optim.minimizer(res)[2])
        u_3 = 1.
        try
        u_3 = find_zero(log_likelihood_ratio_f3, (f_on, 1.))
        catch err
        end
        u_3_fP(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], u_3, rel_div_on)
        res = Optim.optimize(u_3_fP, [m, mu_het])
        m_on[6] = Optim.minimizer(res)[1]*Optim.minimizer(res)[2]*(1-u_3)/u_3
        mu_inc[6] = Optim.minimizer(res)[2]*(1-u_3)/u_3
        del_mu[6] = (1-u_3)*(1+Optim.minimizer(res)[2])
        return [l_1 u_1; l_2 u_2; l_3 u_3; minimum(m_on) maximum(m_on); minimum(mu_inc) maximum(mu_inc); minimum(del_mu) maximum(del_mu)]
    end
end
function CI(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio, m, mu_het, ML)
    function log_likelihood_ratio_1(para)
        if para == m
            return -0.5*chisq_95
        else
            log_likelihood_1(P) = -log_likelihood(mc_p, mc_s, N_ratio, para, P[1])
            res = Optim.optimize(log_likelihood_1, [mu_het])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood(mc_p, mc_s, N_ratio, para, P1) - ML - 0.5*chisq_95
        end
    end
    l_1 = find_zero(log_likelihood_ratio_1, (0., m))
    u_1 = m
    try
        u_1 = find_zero(log_likelihood_ratio_1, (m, maximum([mc_p; mc_s])))
    catch err
        u_1 = find_zero(log_likelihood_ratio_1, (m, 10*maximum([mc_p; mc_s])))
    end
    function log_likelihood_ratio_2(para) 
        if para == mu_het
            return -0.5*chisq_95
        else
            log_likelihood_2(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], para)
            res = Optim.optimize(log_likelihood_2, [m])
            P2 = Optim.minimizer(res)[1]
            return -log_likelihood(mc_p, mc_s, N_ratio, P2, para) - ML - 0.5*chisq_95
        end
    end
    l_2 = 0.
    try 
        l_2 = find_zero(log_likelihood_ratio_2, (0., mu_het))
    catch err
    end
    u_2 = find_zero(log_likelihood_ratio_2, (mu_het, Inf))
    return [l_1 u_1; l_2 u_2]
end
function CI(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio, m, mu_het, f_on, rel_div_on, infer_r::Bool, infer_f::Bool, ML)
    m_on = Vector{Float64}(undef, 8)
    mu_inc = Vector{Float64}(undef, 8)
    del_mu = Vector{Float64}(undef, 8)
    function log_likelihood_ratio_1(para)
        if para == m
            return -0.5*chisq_95
        else
            log_likelihood_1(P) = -log_likelihood(mc_p, mc_s, N_ratio, para, P[1], P[2], P[3])
            res = Optim.optimize(log_likelihood_1, [mu_het, f_on, rel_div_on])
            P1, P2, P3 = Optim.minimizer(res)
            return -log_likelihood(mc_p, mc_s, N_ratio, para, P1, P2, P3) - ML - 0.5*chisq_95
        end
    end
    l_1 = find_zero(log_likelihood_ratio_1, (0., m))
    l_1_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, l_1, P[1], P[2], P[3])
    res = Optim.optimize(l_1_P, [mu_het, f_on, rel_div_on])
    m_on[1] = Optim.minimizer(res)[1]*l_1*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
    mu_inc[1] = Optim.minimizer(res)[1]*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
    del_mu[1] = (1-Optim.minimizer(res)[2])*(1+Optim.minimizer(res)[1])
    u_1 = m
    try
        u_1 = find_zero(log_likelihood_ratio_1, (m, maximum([mc_p; mc_s])))
    catch err
        u_1 = find_zero(log_likelihood_ratio_1, (m, 10*maximum([mc_p; mc_s])))
    end
    u_1_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, u_1, P[1], P[2], P[3])
    res = Optim.optimize(u_1_P, [mu_het, f_on, rel_div_on])
    m_on[2] = Optim.minimizer(res)[1]*u_1*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
    mu_inc[2] = Optim.minimizer(res)[1]*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
    del_mu[2] = (1-Optim.minimizer(res)[2])*(1+Optim.minimizer(res)[1])
    function log_likelihood_ratio_2(para) 
        if para == mu_het
            return -0.5*chisq_95
        else
            log_likelihood_2(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], para, P[2], P[3])
            res = Optim.optimize(log_likelihood_2, [m, f_on, rel_div_on])
            P1, P2, P3 = Optim.minimizer(res)
            return -log_likelihood(mc_p, mc_s, N_ratio, P1, para, P2, P3) - ML - 0.5*chisq_95
        end
    end
    l_2 = 0.
    try 
        l_2 = find_zero(log_likelihood_ratio_2, (0., mu_het))
    catch err
    end
    l_2_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], l_2, P[2], P[3])
    res = Optim.optimize(l_2_P, [m, f_on, rel_div_on])
    m_on[3] = l_2*Optim.minimizer(res)[1]*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
    mu_inc[3] = l_2*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
    del_mu[3] = (1-Optim.minimizer(res)[2])*(1+l_2)
    u_2 = find_zero(log_likelihood_ratio_2, (mu_het, Inf))
    u_2_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], u_2, P[2], P[3])
    res = Optim.optimize(u_2_P, [m, f_on, rel_div_on])
    m_on[4] = u_2*Optim.minimizer(res)[1]*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
    mu_inc[4] = u_2*(1-Optim.minimizer(res)[2])/Optim.minimizer(res)[2]
    del_mu[4] = (1-Optim.minimizer(res)[2])*(1+u_2)
    function log_likelihood_ratio_3(para)
        if para == rel_div_on
            return -0.5*chisq_95
        else 
            log_likelihood_3(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], para, P[3])
            res = Optim.optimize(log_likelihood_3, [m, mu_het, rel_div_on])
            P1, P2, P3 = Optim.minimizer(res)
            return -log_likelihood(mc_p, mc_s, N_ratio, P1, P2, para, P3) - ML - 0.5*chisq_95
        end
    end
    l_3 = 0.
    try
        l_3 = find_zero(log_likelihood_ratio_3, (0., f_on))
    catch err
    end
    l_3_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], l_3, P[3])
    res = Optim.optimize(l_3_P, [m, mu_het, rel_div_on])
    m_on[5] = Optim.minimizer(res)[1]*Optim.minimizer(res)[2]*(1-l_3)/l_3
    mu_inc[5] = Optim.minimizer(res)[2]*(1-l_3)/l_3
    del_mu[5] = (1-l_3)*(1+Optim.minimizer(res)[2])
    u_3 = 1.
    try
        u_3 = find_zero(log_likelihood_ratio_3, (f_on, 1.))
    catch err
    end
    u_3_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], u_3, P[3])
    res = Optim.optimize(u_3_P, [m, mu_het, rel_div_on])
    m_on[6] = Optim.minimizer(res)[1]*Optim.minimizer(res)[2]*(1-u_3)/u_3
    mu_inc[6] = Optim.minimizer(res)[2]*(1-u_3)/u_3
    del_mu[6] = (1-u_3)*(1+Optim.minimizer(res)[2])
    function log_likelihood_ratio_4(para)
        if para == f_on
            return -0.5*chisq_95
        else
            log_likelihood_4(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], P[3], para)
            res = Optim.optimize(log_likelihood_4, [m, mu_het, f_on])
            P1, P2, P3 = Optim.minimizer(res)
            return - log_likelihood(mc_p, mc_s, N_ratio, P1, P2, P3, para) - ML - 0.5*chisq_95
        end
    end
    l_4 = 0.
    try
        l_4 = find_zero(log_likelihood_ratio_4, (0., rel_div_on))
    catch err
    end
    l_4_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], P[3], l_4)
    res = Optim.optimize(l_4_P, [m, mu_het, f_on])
    m_on[7] = Optim.minimizer(res)[1]*Optim.minimizer(res)[2]*(1-Optim.minimizer(res)[3])/Optim.minimizer(res)[3]
    mu_inc[7] = Optim.minimizer(res)[2]*(1-Optim.minimizer(res)[3])/Optim.minimizer(res)[3]
    del_mu[7] = (1-Optim.minimizer(res)[3])*(1+Optim.minimizer(res)[2])
    u_4 = find_zero(log_likelihood_ratio_4, (rel_div_on, Inf))
    u_4_P(P) = -log_likelihood(mc_p, mc_s, N_ratio, P[1], P[2], P[3], u_4)
    res = Optim.optimize(u_4_P, [m, mu_het, f_on])
    m_on[8] = Optim.minimizer(res)[1]*Optim.minimizer(res)[2]*(1-Optim.minimizer(res)[3])/Optim.minimizer(res)[3]
    mu_inc[8] = Optim.minimizer(res)[2]*(1-Optim.minimizer(res)[3])/Optim.minimizer(res)[3]
    del_mu[8] = (1-Optim.minimizer(res)[3])*(1+Optim.minimizer(res)[2])
    return [l_1 u_1; l_2 u_2; l_3 u_3; l_4 u_4; minimum(m_on) maximum(m_on); minimum(mu_inc) maximum(mu_inc); minimum(del_mu) maximum(del_mu)]
end

function CV(mc::Vector{Int}, fit_m)
	mc_9010 = copy(mc)
    scores = zeros(Float64, 10)
    n = Int(length(mc)/10)
    for i = 1:10
        c = mc_9010[1:n]
        deleteat!(mc_9010, 1:n) 
        log_likelihood_para(para) = -log_likelihood(mc_9010, para[1], fit_m)
	    res = Optim.optimize(log_likelihood_para, [initial_m(mc_9010, 1000)])
	    if Optim.converged(res) == true
            scores[i] = -log_likelihood(c, Optim.minimizer(res)[1], fit_m)                                                                            
	    end 
        mc_9010 = [mc_9010; c]
	end
    return scores
end 
function CV(mc::Vector{Int})
    mc_9010 = copy(mc)
	scores = zeros(Float64, 10)
    n = Int(length(mc)/10)
    for i = 1:10
        c = mc_9010[1:n]
        deleteat!(mc_9010, 1:n)                                                     
	    log_likelihood_para(para) = -log_likelihood(mc_9010, para[1], para[2])              
	    res = Optim.optimize(log_likelihood_para, [initial_m(mc_9010, 1000), 1.]) 
	    if Optim.converged(res) == true
		    p = Optim.minimizer(res)
            scores[i] = -log_likelihood(c, p[1], p[2])                         
        end
        mc_9010 = [mc_9010; c]
	end
	return scores                                                        
end
function CV(mc_p::Vector{Int}, mc_s::Vector{Int})
	mc_p_9010 = copy(mc_p)
    mc_s_9010 = copy(mc_s)
	scores = zeros(Float64, 10)
    n_p = Int(length(mc_p)/10)
    n_s = Int(length(mc_s)/10)
    for i = 1:10
        c_p = mc_p_9010[1:n_p]
        c_s = mc_s_9010[1:n_s]
        deleteat!(mc_p_9010, 1:n_p)        
        deleteat!(mc_s_9010, 1:n_s)                                                          
	    log_likelihood_para(para) = -log_likelihood(mc_p_9010, para[1], mc_s_9010, para[2], para[3])
	    res = Optim.optimize(log_likelihood_para, [1., 1., 1.])                                      
	    if Optim.converged(res) == true
		    p = Optim.minimizer(res)
            scores[i] = -log_likelihood(c_p, p[1], c_s, p[2], p[3])
        end    
        mc_p_9010 = [mc_p_9010; c_p]
        mc_s_9010 = [mc_s_9010; c_s]
    end
	return scores
end
function CV(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio, f_on, rel_div_on)
	mc_p_9010 = copy(mc_p)
    mc_s_9010 = copy(mc_s)
	scores = zeros(Float64, 10)
    n_p = Int(length(mc_p)/10)
    n_s = Int(length(mc_s)/10)
    for i = 1:10
        c_p = mc_p_9010[1:n_p]
        c_s = mc_s_9010[1:n_s]
        deleteat!(mc_p_9010, 1:n_p)        
        deleteat!(mc_s_9010, 1:n_s)
        log_likelihood_para(para) = -log_likelihood(mc_p_9010, mc_s_9010, N_ratio, para[1], para[2], f_on, rel_div_on)
        res = Optim.optimize(log_likelihood_para, [1., 1.])                                                 
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            scores[i] = -log_likelihood(c_p, c_s, N_ratio, p[1], p[2], f_on, rel_div_on)
        end
        mc_p_9010 = [mc_p_9010; c_p]
        mc_s_9010 = [mc_s_9010; c_s]
    end
    return scores
end
function CV(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio, f_on, rel_div_on, infer_r_f::Bool)
    mc_p_9010 = copy(mc_p)
    mc_s_9010 = copy(mc_s)
	scores = zeros(Float64, 10)
    n_p = Int(length(mc_p)/10)
    n_s = Int(length(mc_s)/10)
    if infer_r_f == true
        for i = 1:10
            c_p = mc_p_9010[1:n_p]
            c_s = mc_s_9010[1:n_s]
            deleteat!(mc_p_9010, 1:n_p)        
            deleteat!(mc_s_9010, 1:n_s) 
            log_likelihood_para_r(para) = -log_likelihood(mc_p_9010, mc_s_9010, N_ratio, para[1], para[2], f_on, para[3])
            res = Optim.optimize(log_likelihood_para_rp, [1., 1., rel_div_on])                                                 
            if Optim.converged(res) == true
                p = Optim.minimizer(res)
                scores[i] = -log_likelihood(c_p, c_s, N_ratio, p[1], p[2], f_on, p[3])
            end
            mc_p_9010 = [mc_p_9010; c_p]
            mc_s_9010 = [mc_s_9010; c_s]
        end
    else
        for i = 1:10
            c_p = mc_p_9010[1:n_p]
            c_s = mc_s_9010[1:n_s]
            deleteat!(mc_p_9010, 1:n_p)        
            deleteat!(mc_s_9010, 1:n_s) 
            log_likelihood_para_f(para) = -log_likelihood(mc_p_9010, mc_s_9010, N_ratio, para[1], para[2], para[3], rel_div_on)
            res = Optim.optimize(log_likelihood_para_f, [1., 1., f_on])                                                 
            if Optim.converged(res) == true
                p = Optim.minimizer(res)
                scores[i] = -log_likelihood(c_p, c_s, N_ratio, p[1], p[2], p[3], rel_div_on)
            end
            mc_p_9010 = [mc_p_9010; c_p]
            mc_s_9010 = [mc_s_9010; c_s]
        end
    end
    return scores
end
function CV(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio)
    mc_p_9010 = copy(mc_p)
    mc_s_9010 = copy(mc_s)
    scores = zeros(Float64, 10)
    n_p = Int(length(mc_p)/10)
    n_s = Int(length(mc_s)/10)
    for i = 1:10
        c_p = mc_p_9010[1:n_p]
        c_s = mc_s_9010[1:n_s]
        deleteat!(mc_p_9010, 1:n_p)        
        deleteat!(mc_s_9010, 1:n_s)     
        log_likelihood_para(para) = -log_likelihood(mc_p_9010, mc_s_9010, N_ratio, para[1], para[2])     
        res = Optim.optimize(log_likelihood_para, [1., 1.])                     
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            scores[i] = -log_likelihood(c_p, c_s, N_ratio, p[1], p[2]) 
        end
        mc_p_9010 = [mc_p_9010; c_p]
        mc_s_9010 = [mc_s_9010; c_s]
    end
    return scores
end
function CV(mc_p::Vector{Int}, mc_s::Vector{Int}, N_ratio, f_on, infer_r_f::Bool)
    mc_p_9010 = copy(mc_p)
    mc_s_9010 = copy(mc_s)
    scores = zeros(Float64, 10)
    n_p = Int(length(mc_p)/10)
    n_s = Int(length(mc_s)/10)
    for i = 1:10
        c_p = mc_p_9010[1:n_p]
        c_s = mc_s_9010[1:n_s]
        deleteat!(mc_p_9010, 1:n_p)        
        deleteat!(mc_s_9010, 1:n_s)  
        log_likelihood_para(para) = -log_likelihood(mc_p_9010, mc_s_9010, N_ratio, para[1], para[2], para[3], para[4])
        res = Optim.optimize(log_likelihood_para, [1., 1., f_on, 0.])                                                 
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            scores[i] = -log_likelihood(c_p, c_s, N_ratio, p[1], p[2], p[3], p[4])
        end
        mc_p_9010 = [mc_p_9010; c_p]
        mc_s_9010 = [mc_s_9010; c_s]
    end
    return scores
end

# Probability generating function method used to set initial values of the maximum likelihood estimation, based on
# Gillet-Markowska, A., Louvel, G., & Fischer, G. (2015). bz-rates: A web tool to estimate mutation rates from fluctuation analysis. G3: Genes, Genomes, Genetics, 5(11), 2323–2327. https://doi.org/10.1534/g3.115.019836
function empirical_pgf(z, x) # Empirical probability generating function calculated from observed data
    g = 0
    for i in x
        g += z^i
    end
    g /= length(x)
    return g
end
function initial_m(z, mc)                     # Estimate number of mutations for a homogeneous population and given z        
    if z == 0.
        return log(empirical_pgf(z, mc)) 
    else
        return z/((1-z)*log(1-z)) * log(empirical_pgf(z, mc))
    end
end
function initial_m(mc, z_values::Int)         # Estimate number of mutations for a homogeneous population by averaging over a number of z values
    m = 0.
    for i = 0:z_values-1
        m += initial_m(i/z_values, mc)
    end
    return maximum([m/z_values, 0.])
end
function initial_mu_het(z, mc, m)             # Estimate the mutation-rate heterogeneity for the heterogeneous-response model for given z   
    if z == 0.
        return -log(empirical_pgf(z, mc))/m - 1
    else
        return -log(empirical_pgf(z, mc))/(m*(1-z)) + log(1-z)/z
    end
end
function initial_mu_het(mc, m, z_values::Int) # Estimate the mutation-rate heterogeneity for the heterogeneous-response model by averaging over a number of z values 
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