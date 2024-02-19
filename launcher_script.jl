import Pkg
#Pkg.activate("packages")
#Pkg.instantiate()
#include("population_dynamics.jl")
#include("inference.jl")
try
    mkdir("output_data")
catch e
end
try
    mkdir("inferred_parameters")
catch e
end
using CSV, DataFrames, Random

# Simulate fluctuation assays for a parameter range given in the file of input parameters p.csv
function simulate_fluctuation_assays(range, range_2=""; set_seed=false, S1=false)
    N0 = 10^4
    expected_M = 1
    division_off = 1.
    fitness_m_p, fitness_m_s = 1., 1.
    mutation_off_p, mutation_off_s = 10^-8., 10^-8.
    death_off = 0.
    division_on = 0.
    death_on = 0.
    Nf = 10^9
    num_cultures = 10^4
    p = DataFrame(CSV.File("input_parameters/"*range*".csv"))   # Reading other input parameters             
    f0_on = p.f0_on[1]
    f_stat = false
    switching = p.alpha[1]
    mutation_on = p.nu_on[1]
    r_parameter = p.r_parameter[1]                              # This parameter is varied; the other parameters are fixed
    r_start = p.r_start[1]
    r_end = p.r_end[1]
    r_i = p.r_increment[1]
    if set_seed == true
        rng=Random.seed!(p.random_seed[1])                      # Seed for random number generator to reproduce the exact data in the manuscript
    end
    if S1 == true
        n = Matrix{Int}(undef, (100,2))
        n_t1 = DataFrame()
        n_tf = DataFrame()
    else
        mutant_counts = DataFrame()
        p_final = DataFrame()
    end
    if range_2 == ""                                            # Only one parameter is varied      
        r2_parameter = ""
        r2_start = 1
        r2_i = 1
        r2_end = 1         
    else                                                        # A second parameter is varied making the parameter range 2D
        p2 = DataFrame(CSV.File("input_parameters/"*range_2*".csv"))    
        r2_parameter = p2.r_parameter[1]                              
        r2_start = p2.r_start[1]
        r2_end = p2.r_end[1]
        r2_i = p2.r_increment[1]
        if p2.f0_on[1] == "stat"
            f_stat = true
        end
        j2 = 1  
        try
            mkdir("output_data/"*range_2)
        catch e
        end
        range_2 *= "/"
    end    
    for r2 = r2_start:r2_i:r2_end
        if r2_parameter == "log_alpha"
            switching = 10^r2 
        elseif r2_parameter == "delta_on"
            death_on = r2
        elseif r2_parameter == "rho"
            fitness_m_p, fitness_m_s = r2, r2
        elseif r2_parameter == "rho_s"
            fitness_m_s = r2
        end
        if f_stat
            f0_on = switching/(division_off-death_off)
        end
        j = 1
        for r = r_start:r_i:r_end
            if r_parameter == "log_nu_on"
                mutation_on = 10^r
            elseif r_parameter == "delta_off"
                death_off = r
            elseif r_parameter == "gamma_on"
                division_on = r
            elseif r_parameter == "log_nu_off_s"
                mutation_off_s = 10^r
            elseif r_parameter == "rho"
                fitness_m_p, fitness_m_s = r, r
            end
            if S1 == true
                t1 = t_first_m(N0, division_off-switching-death_off, mutation_off_s*(1-f0_on)+mutation_on*f0_on)
                tf = t_final(N0, division_off-switching-death_off, Nf)
                for i = 1:size(n)[1]
                    n[i,1] = non_mutants_on(t1, N0*(1-f0_on), division_off-switching-death_off, switching, Int(round(N0*f0_on)), division_on, death_on, N0)
                    n[i,2] = non_mutants_on(tf-t1, pop_size(t1, N0*(1-f0_on), division_off-switching-death_off), division_off-switching-death_off, switching, n[i,1], division_on, death_on, N0)
                end
                n_t1[:,"$j"] = n[:,1]
                n_tf[:,"$j"] = n[:,2]
            else
                T = t_expected_m(N0*(1-f0_on), division_off-switching-death_off, mutation_off_s, switching, N0*f0_on, division_on, mutation_on, expected_M)
                mc, Nf, f_on = mutant_count(T, N0, division_off, mutation_off_s, switching, division_on, mutation_on, num_cultures, fitness_m_off=fitness_m, death_off=death_off, death_on=death_on, f0_on=f0_on)
                mutant_counts[:,"$j"] = mc
                p_final[:,"$j"] = [Nf, f_on, T]
            end
            j += 1
        end
        if S1 == true
            CSV.write("output_data/"*range_2*"n_t1-"*range*"_$(r2_parameter)_$j2.csv", n_t1)
            CSV.write("output_data/"*range_2*"n_tf-"*range*"_$(r2_parameter)_$j2.csv", n_tf)
            j2 += 1
        else
            T_p = t_expected_m(N0, division_off, mutation_off_p, 0., 0, 0., 0., expected_M)
            mc_p, Nf_p = mutant_count(T_p, N0, division_off, mutation_off_p, num_cultures, fitness_m=fitness_m)
            mutant_counts.p = mc_p
            p_final.p = [Nf_p, 0., T_p]
            if range_2 == "" 
                CSV.write("output_data/mutant_counts-"*range*".csv", mutant_counts)
                CSV.write("output_data/p_final-"*range*".csv", p_final)
            else
                CSV.write("output_data/"*range_2*"mutant_counts-"*range*"_$(r2_parameter)_$j2.csv", mutant_counts)
                CSV.write("output_data/"*range_2*"p_final-"*range*"_$(r2_parameter)_$j2.csv", p_final)
                j2 += 1
            end
        end
    end
end

# Mutation rate estimation for ranges of parameters
# The inferred parameters is saved into a folder named p (more information can be found in the README file)
function infer_mutation_rates(p, m; p_folder="") # Parameter range p and inference method m; optional input when a second parameter range is considered
    parameters = DataFrame(CSV.File("input_parameters/"*p*".csv"))
    r_parameter = parameters.r_parameter[1]                # The parameter that is varied 
    R = parameters.number_fluctuation_assays[1]
    mc_bound = 1000                                                # Discarding all mutant counts >mc_bound
    try
        mkdir("inferred_parameters/"*p)                            # Make a folder to store the inferred parameters in
    catch e
    end
    p = p*p_folder                                                 # In the case of a second parameter range, a subfolder needs to be made (see function below)
    try
        mkdir("inferred_parameters/"*p)                            # Make a folder to store the inferred parameters in
    catch e
    end
    mc_p = DataFrame(CSV.File("output_data/"*p*"/mutant_counts_p.csv"))             # Read the output data without stress
    J = Int(size(mc_p)[2]/R)
    Nf_p = DataFrame(CSV.File("output_data/"*p*"/population_sizes_p.csv")).value[1]
    # For some methods the relative division rate of response-on cells is set to the true value
    if r_parameter == "divisions_on"
        rel_div_on = collect(parameters.r_start[1]:parameters.r_increment[1]:parameters.r_end[1]) ./ parameters.divisions_off
    else
        rel_div_on = fill(parameters.divisions_on[1], J)
    end 
    for j = 1:J
        mc_s = DataFrame(CSV.File("output_data/"*p*"/mutant_counts_$j.csv"))    # Read the output data with stress
        pop_sizes = DataFrame(CSV.File("output_data/"*p*"/population_sizes_$j.csv")) 
        Nf_s = pop_sizes.value[1]
        f_on = pop_sizes.value[2]
        inferred_para = DataFrame()
        # Inferred parameters are listed below. Parameters that are not inferred by a particular method are simply left to be zero
        mu_offs = zeros(Float64, R)
        mu_ons = zeros(Float64, R)
        div_ons = zeros(Float64, R)
        f_ons = zeros(Float64, R)
        mu_p = zeros(Float64, R)
        fitm_p = zeros(Float64, R)
        mu_s = zeros(Float64, R)
        fitm_s = zeros(Float64, R)
        AICs = zeros(Float64, R)
        # Inference methods are listed below
        # Inference under the heterogeneous-response model for known fraction of response-on subpopulation and with the relative division rate of response-on cells as an inference parameter 
        if m == "het_infer_div" 
            for i = 1:R
                mu_offs[i], mu_ons[i], div_ons[i], AICs[i] = estimate_mu_het(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, f_on, rel_div_on="infer")
            end
        # Inference under the heterogeneous-response model for known fraction of response-on subpopulation and setting the relative division rate of response-on cells to the true value
        elseif m == "het_set_div"
            for i = 1:R
                mu_offs[i], mu_ons[i], div_ons[i], AICs[i] = estimate_mu_het(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, f_on, rel_div_on=rel_div_on[j])
            end
        # Inference under the heterogeneous-response model for known fraction of response-on subpopulation and setting the relative division rate of response-on cells to zero
        elseif m == "het_zero_div"
            for i = 1:R
                mu_offs[i], mu_ons[i], div_ons[i], AICs[i] = estimate_mu_het(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, f_on)
            end
        # Inference under the heterogeneous-response model for unknown fraction of response-on subpopulation
        elseif m == "het_unknown_fraction"
            for i = 1:R
                mu_offs[i], mu_ons[i], f_ons[i], div_ons[i], AICs[i] = estimate_mu_het(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s)
            end
        # Inference under the heterogeneous-response model for unknown fraction of response-on subpopulation and with the relative division rate of response-on cells as an inference parameter
        elseif m == "het_unknown_fraction_infer_div"
            for i = 1:R
                mu_offs[i], mu_ons[i], f_ons[i], div_ons[i], AICs[i] = estimate_mu_het(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, rel_div_on="infer")
            end
        # Inference under homogeneous-response model with the differential fitness of mutants as an inference parameter
        elseif m == "hom_infer_fit"
            for i = 1:R
                mu_p[i], fitm_p[i], mu_s[i], fitm_s[i], AICs[i] = estimate_mu_hom(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, fit_m="infer")
            end
        # Inference under homogeneous-response model with the differential fitness of mutants as a joint inference parameter
        elseif m == "hom_joint_fit"
            for i = 1:R
                mu_p[i], fitm_p[i], mu_s[i], fitm_s[i], AICs[i] = estimate_mu_hom(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, fit_m="joint")
            end
        # Inference under homogeneous-response model without the differential fitness of mutants as an inference parameter (set to one)
        elseif m == "hom_no_fit"
            for i = 1:R
                mu_p[i], fitm_p[i], mu_s[i], fitm_s[i], AICs[i] = estimate_mu_hom(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s)
            end
        end
        inferred_para[:, "mutation_rate_off"] = mu_offs
        inferred_para[:, "mutation_rate_on"] = mu_ons
        inferred_para[:, "division_rate_on"] = div_ons
        inferred_para[:, "fraction_on"] = f_ons
        inferred_para[:, "mutation_rate_p"] = mu_p
        inferred_para[:, "fitness_mants_p"] = fitm_p
        inferred_para[:, "mutation_rate_s"] = mu_s
        inferred_para[:, "fitness_mants_s"] = fitm_s
        inferred_para[:, "AIC"] = AICs
        CSV.write("inferred_parameters/"*p*"/"*m*"_$j.csv", inferred_para) # Save data frame with inferred parameters as csv file
    end
end

function infer_mutation_rates(p, m, p2) # Two parameter ranges are given as input
    parameters2 = DataFrame(CSV.File("input_parameters/"*p2*".csv"))                            # Read the input parameters of the second parameter range
    r_parameter2 = parameters2.r_parameter[1]                                           # The second parameter that is varied 
    J2 = Int(round((parameters2.r_end[1]-parameters2.r_start[1])/parameters2.r_increment[1]))+1 # in J2 steps
    # For each step, a subfolder will be made with the name listed below
    if r_parameter2 == "mutation_increase_log"
        p_folder = "/mu-inc_"   
    elseif r_parameter2 == "divisions_on"
        p_folder = "/rel-div_" 
    elseif r_parameter2 == "fitness_mants"
        p_folder = "/fit-mutants_" 
    elseif r_parameter2 == "deaths_off"
        p_folder = "/death-off_" 
    elseif r_parameter2 == "deaths_on"
        p_folder = "/death-on_" 
    elseif r_parameter2 == "switchings_log"
        p_folder = "/switch_" 
    end
    for j2 = 1:J2
        infer_mutation_rates(p, m, p_folder=p_folder*"$j2") 
    end
end

# Reproduce all data in the manuscript

function data_inference_manuscript()
    # Parameter regime: mutation-rate increase x relative switching rate
    # Estimation method: heterogeneous-response model with setting the relative division rate of response-on cells to zero (known fraction of response-on subpopulation)
    simulate_fluctuation_assays("range-nu_on", "range-alpha-f0", set_seed=true)
    infer_mutation_rates("range_mu-inc", "het_zero_div", "range_switching")

    # Parameter regime: death rate of response-off x -on cells, for switching rates 0.01 and 0.05
    # Estimation method: heterogeneous-response model with setting the relative division rate of response-on cells to zero (known fraction of response-on subpopulation)
    for i in [1, 5]
        for j in ["_off", "_on", ""]
            simulate_fluctuation_assays("range-delta"*j*"-alpha$i", set_seed=true)
            infer_mutation_rates("range_death-off_switch-$i", "het_zero_div", "range_death-on_switch-$i")
        end
    end

    # Parameter regime: differential fitness of response-off mutants
    # Estimation method: heterogeneous-response model with setting the relative division rate of response-on cells to zero (known fraction of response-on subpopulation)
    simulate_fluctuation_assays("range-mu_off", "range-rho", set_seed=true)
    infer_mutation_rates("range_fit-mut", "het_zero_div")

    # Parameter regime: relative division rate of response-on cells
    # Estimation methods
    # (i) Heterogeneous-response model with setting the relative division rate of response-on cells to zero/true value or inferring it (known fraction of response-on subpopulation)
    # (ii) Heterogeneous-response model with setting the relative division rate of response-on cells to zero or inferring it (unknown fraction of response-on subpopulation)
    # (iii) Homogeneous-response model without/with/jointly inferring the differential fitness of mutants
    for i in [10, 100]
        simulate_fluctuation_assays("range-gamma_on-increase$i", set_seed=true)
        for m in ["het_zero_div", "het_set_div", "het_infer_div", "het_unknown_fraction", "het_unknown_fraction_infer_div", "hom_no_fit", "hom_infer_fit", "hom_joint_fit"]
            infer_mutation_rates(p, m)
        end
    end
end

function data_supplementary_material()
    # Parameter regime: switching rate x relative division rate of response-on cells
    # Simulating response-on non-mutants stochastically to test assumption S1
    for j in ["f0", "fstat"]
        simulate_fluctuation_assays("range-gamma_on-increase100", "range-alpha-"*j, set_seed=true, S1=true)
    end
end
