import Pkg
Pkg.activate("packages")
include("population_dynamics.jl")
include("inference.jl")
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
# The output data is saved into a folder named p (more information can be found in the README file)
function simulate_fluctuation_assays(p; p2="", set_seed=false, A1=false)
    # Reading the input parameters
    parameters = DataFrame(CSV.File("input_parameters/"*p*".csv"))             
    if set_seed == true
        rng=Random.seed!(parameters.random_seed[1])        # Seed for random number generator to reproduce the exact data in the manuscript
    end
    try
        mkdir("output_data/"*p)                            # Make a folder to store the output data in
    catch e
    end
    N0 = parameters.initial_population_size[1]
    division_off = parameters.divisions_off[1]
    mutation_off = parameters.mutations_off[1]
    switching = parameters.switchings[1]
    division_on = parameters.divisions_on[1]
    num_cultures = parameters.number_cultures[1]
    mutation_on = parameters.mutations_on[1]
    expected_M = parameters.expected_mutations[1]
    R = parameters.number_fluctuation_assays[1]                                # Number of fluctuation assays that is simulated
    range_parameter = parameters.range_parameter[1]                            # This parameter is varied; the other parameters are fixed
    r_start = parameters.r_start[1]
    r_end = parameters.r_end[1]
    r_i = parameters.r_increment[1]
    j = 1
    fitness_m_off = 1.
    death_off = 0.
    death_on = 0.
    f0_on = switching/(division_off-death_off)
    Nf = 10^9
    # Data frames in which the simulated data will be stored
    if A1 == true
        f0_on = [0., 0.01, switching/(division_off-death_off)]
        f = length(f0_on)
        popsize_nm_on_t1 = DataFrame()
        popsize_nm_on_tf = DataFrame()
        nm1 = Matrix{Int}(undef, (R, f))
        nmf = Matrix{Int}(undef, (R, f))
    else
        mutant_counts = DataFrame()
        population_sizes = DataFrame(parameter=["final_pop_size", "fraction_subpop_on"]) # Under stress, final population size and fraction of response-on subpopulation are stored
        mutant_counts_p = DataFrame()
        population_sizes_p = DataFrame(parameter=["final_pop_size"])
        T_p = t_expected_m(N0, division_off, mutation_off, 0., 0, 0., 0., expected_M) # Under no stress (ns) conditions, there is no response-on subpopulation
        Nf_p = pop_size(T_p, N0, division_off)
        population_sizes_p.value = [Nf_p]
    end
    # Default: only one parameter is varied
    if p2 == ""                                                         
        for r = r_start:r_i:r_end
            if range_parameter == "mutation_increase_log"
                mutation_on = 10^r * mutation_off   
            elseif range_parameter == "divisions_on"
                division_on = r
            elseif range_parameter == "fitness_mutants"
                fitness_m_off = r
            elseif range_parameter == "deaths_off"
                death_off = r
            elseif range_parameter == "deaths_on"
                death_on = r
            elseif range_parameter == "switchings_log"
                switching = 10^r
            end
            if A1 == true
                f0_on[end] = switching/(division_off-death_off)
                tf = t_final(N0, division_off-switching-death_off, Nf)
                for f = 1:length(f0_on)
                    t1 = t_first_m(N0, division_off-switching-death_off, mutation_off*(1-f0_on[f])+mutation_on*f0_on[f])
                    for i = 1:R
                        nm1[i, f] = non_mutants_on(t1, N0*(1-f0_on[f]), division_off-switching-death_off, switching, Int(round(N0*f0_on[f])), division_on, death_on, N0)
                        nmf[i, f] = non_mutants_on(tf-t1, pop_size(t1, N0*(1-f0_on[f]), division_off-switching-death_off), division_off-switching-death_off, switching, nm1[i, f], division_on, death_on, N0)
                    end
                    popsize_nm_on_t1[:, "$f"] = nm1[: ,f]
                    popsize_nm_on_tf[:, "$f"] = nmf[: ,f]
                end
                CSV.write("output_data/"*p*"/popsize_nm_on_t1_$j.csv", popsize_nm_on_t1)
                CSV.write("output_data/"*p*"/popsize_nm_on_tf_$j.csv", popsize_nm_on_tf)
            else
                T = t_expected_m(N0*(1-f0_on), division_off-switching-death_off, mutation_off, switching, N0*f0_on, division_on, mutation_on, expected_M)
                for i = 1:R
                    mc, Nf, f_on = mutant_count(T, N0, division_off, mutation_off, switching, division_on, mutation_on, num_cultures, fitness_m_off=fitness_m_off, death_off=death_off, death_on=death_on, f0_on=f0_on)
                    mutant_counts[:, "$i"] = mc
                    mc_p, Nf_p = mutant_count(T_p, N0, division_off, mutation_off, num_cultures)
                    mutant_counts_p[:, "$(R*(j-1)+i)"] = mc_p
                end
                mc, Nf, f_on = mc, Nf, f_on = mutant_count(T, N0, division_off, mutation_off, switching, division_on, mutation_on, num_cultures, fitness_m_off=fitness_m_off, death_off=death_off, death_on=death_on, f0_on=f0_on)
                population_sizes.value = [Nf, f_on]
                CSV.write("output_data/"*p*"/mutant_counts_$j.csv", mutant_counts)
                CSV.write("output_data/"*p*"/population_sizes_$j.csv", population_sizes)
            end
            j += 1
        end
        if A1 == true
        else
            CSV.write("output_data/"*p*"/mutant_counts_p.csv", mutant_counts_p)
            CSV.write("output_data/"*p*"/population_sizes_p.csv", population_sizes_p)
        end
    # A second parameter is varied making the parameter range 2D
    else                                                                   
        parameters2 = DataFrame(CSV.File("input_parameters/"*p2*".csv"))
        range_parameter2 = parameters2.range_parameter[1]                            
        r_start2 = parameters2.r_start[1]
        r_end2 = parameters2.r_end[1]
        r_i2 = parameters2.r_increment[1]
        j2 = 1
        for r2 = r_start2:r_i2:r_end2
            if range_parameter2 == "mutation_increase_log"
                mutation_on = 10^r2 * mutation_off
                p_folder = "/mu-inc_"   
            elseif range_parameter2 == "divisions_on"
                division_on = r2
                p_folder = "/rel-div_" 
            elseif range_parameter2 == "fitness_mutants"
                fitness_m_off = r2
                p_folder = "/fit-mutants_" 
            elseif range_parameter2 == "deaths_off"
                death_off = r2
                p_folder = "/death-off_" 
            elseif range_parameter2 == "deaths_on"
                death_on = r2
                p_folder = "/death-on_" 
            elseif range_parameter2 == "switchings_log"
                switching = 10^r2
                p_folder = "/switch_" 
            end
            j = 1
            try
                mkdir("output_data/"*p*p_folder*"$j2")
            catch e
            end
            for r = r_start:r_i:r_end
                if range_parameter == "mutation_increase_log"
                    mutation_on = 10^r * mutation_off   
                elseif range_parameter == "divisions_on"
                    division_on = r
                elseif range_parameter == "fitness_mutants"
                    fitness_m_off = r
                elseif range_parameter == "deaths_off"
                    death_off = r
                elseif range_parameter == "deaths_on"
                    death_on = r
                elseif range_parameter == "switchings_log"
                    switching = 10^r
                end
                if A1 == true
                    f0_on[end] = switching/(division_off-death_off)
                    tf = t_final(N0, division_off-switching-death_off, Nf)
                    for f = 1:length(f0_on)
                        t1 = t_first_m(N0, division_off-switching-death_off, mutation_off*(1-f0_on[f])+mutation_on*f0_on[f])
                        for i = 1:R
                            nm1[i, f] = non_mutants_on(t1, N0*(1-f0_on[f]), division_off-switching-death_off, switching, Int(round(N0*f0_on[f])), division_on, death_on, N0)
                            nmf[i, f] = non_mutants_on(tf-t1, pop_size(t1, N0*(1-f0_on[f]), division_off-switching-death_off), division_off-switching-death_off, switching, nm1[i, f], division_on, death_on, N0)
                        end
                        popsize_nm_on_t1[:, "$f"] = nm1[: ,f]
                        popsize_nm_on_tf[:, "$f"] = nmf[: ,f]
                    end
                    CSV.write("output_data/"*p*p_folder*"$j2/popsize_nm_on_t1_$j.csv", popsize_nm_on_t1)
                    CSV.write("output_data/"*p*p_folder*"$j2/popsize_nm_on_tf_$j.csv", popsize_nm_on_tf)
                else
                    T = t_expected_m(N0*(1-f0_on), division_off-switching-death_off, mutation_off, switching, N0*f0_on, division_on, mutation_on, expected_M)
                    for i = 1:R
                        mc, Nf, f_on = mc, Nf, f_on = mutant_count(T, N0, division_off, mutation_off, switching, division_on, mutation_on, num_cultures, fitness_m_off=fitness_m_off, death_off=death_off, death_on=death_on, f0_on=f0_on)
                        mutant_counts[:, "$i"] = mc
                        mc_p, Nf_p = mutant_count(T_p, N0, division_off, mutation_off, num_cultures)
                        mutant_counts_p[:, "$(R*(j-1)+i)"] = mc_p
                    end
                    mc, Nf, f_on = mc, Nf, f_on = mutant_count(T, N0, division_off, mutation_off, switching, division_on, mutation_on, num_cultures, fitness_m_off=fitness_m_off, death_off=death_off, death_on=death_on, f0_on=f0_on)
                    population_sizes.value = [Nf, f_on]
                    CSV.write("output_data/"*p*p_folder*"$j2/mutant_counts_$j.csv", mutant_counts)
                    CSV.write("output_data/"*p*p_folder*"$j2/population_sizes_$j.csv", population_sizes)
                end
                j += 1
            end
            if A1 == true
            else
                CSV.write("output_data/"*p*p_folder*"$j2/mutant_counts_p.csv", mutant_counts_p)
                CSV.write("output_data/"*p*p_folder*"$j2/population_sizes_p.csv", population_sizes_p)
            end
            j2 += 1
        end
    end
end

# Mutation rate estimation for ranges of parameters
# The inferred parameters is saved into a folder named p (more information can be found in the README file)
function infer_mutation_rates(p, m; p_folder="") # Parameter range p and inference method m; optional input when a second parameter range is considered
    parameters = DataFrame(CSV.File("input_parameters/"*p*".csv"))
    range_parameter = parameters.range_parameter[1]                                         # The parameter that is varied 
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
    if range_parameter == "divisions_on"
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
                mu_offs[i], mu_ons[i], AICs[i] = estimate_mu_het(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, f_on, rel_div_on=rel_div_on[j])
            end
        # Inference under the heterogeneous-response model for known fraction of response-on subpopulation and setting the relative division rate of response-on cells to zero
        elseif m == "het_zero_div"
            for i = 1:R
                mu_offs[i], mu_ons[i], AICs[i] = estimate_mu_het(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, f_on)
            end
        # Inference under the heterogeneous-response model for unknown fraction of response-on subpopulation and with the relative division rate of response-on cells as an inference parameter 
        elseif m == "het_infer_div_unknown_fraction"
            for i = 1:R
                mu_offs[i], mu_ons[i], div_ons[i], f_ons[i], AICs[i] = estimate_mu_het(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, rel_div_on="infer")
            end
        # Inference under the heterogeneous-response model for unknown fraction of response-on subpopulation and setting the relative division rate of response-on cells to zero
        elseif m == "het_zero_div_unknown_fraction"
            for i = 1:R
                mu_offs[i], mu_ons[i], f_ons[i], AICs[i] = estimate_mu_het(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s)
            end
        # Inference under the heterogeneous-response model for unknown fraction of response-on subpopulation and setting the relative division rate of response-on cells to the true value
        elseif m == "het_set_div_unknown_fraction"
            for i = 1:R
                mu_offs[i], mu_ons[i], f_ons[i], AICs[i] = estimate_mu_het(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, rel_div_on=rel_div_on[j])
            end
        # Inference under homogeneous-response model with the differential fitness of mutants as an inference parameter
        elseif m == "hom_infer_fit"
            for i = 1:R
                mu_p[i], fitm_p[i], mu_s[i], fitm_s[i], AICs[i] = estimate_mu_hom(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, fit_m="infer")
            end
        # Inference under homogeneous-response model with the differential fitness of mutants as a joint inference parameter
        elseif m == "hom_joint_fit"
            for i = 1:R
                mu_p[i], mu_s[i], fitm, AICs[i] = estimate_mu_hom(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s, fit_m="joint")
                fitm_s[i] = fitm
                fitm_p[i] = fitm
            end
        # Inference under homogeneous-response model without the differential fitness of mutants as an inference parameter (set to one)
        elseif m == "hom_no_fit"
            for i = 1:R
                mu_p[i], mu_s[i], AICs[i] = estimate_mu_hom(mc_p[:,R*(j-1)+i][mc_p[:,R*(j-1)+i].<mc_bound], Nf_p, mc_s[:,i][mc_s[:,i].<mc_bound], Nf_s)
            end
        end
        inferred_para[:, "mutation_rate_off"] = mu_offs
        inferred_para[:, "mutation_rate_on"] = mu_ons
        inferred_para[:, "division_rate_on"] = div_ons
        inferred_para[:, "fraction_on"] = f_ons
        inferred_para[:, "mutation_rate_p"] = mu_p
        inferred_para[:, "fitness_mutants_p"] = fitm_p
        inferred_para[:, "mutation_rate_s"] = mu_s
        inferred_para[:, "fitness_mutants_s"] = fitm_s
        inferred_para[:, "AIC"] = AICs
        CSV.write("inferred_parameters/"*p*"/"*m*"_$j.csv", inferred_para) # Save data frame with inferred parameters as csv file
    end
end

function infer_mutation_rates(p, m, p2) # Two parameter ranges are given as input
    parameters2 = DataFrame(CSV.File("input_parameters/"*p2*".csv"))                            # Read the input parameters of the second parameter range
    range_parameter2 = parameters2.range_parameter[1]                                           # The second parameter that is varied 
    J2 = Int(round((parameters2.r_end[1]-parameters2.r_start[1])/parameters2.r_increment[1]))+1 # in J2 steps
    # For each step, a subfolder will be made with the name listed below
    if range_parameter2 == "mutation_increase_log"
        p_folder = "/mu-inc_"   
    elseif range_parameter2 == "divisions_on"
        p_folder = "/rel-div_" 
    elseif range_parameter2 == "fitness_mutants"
        p_folder = "/fit-mutants_" 
    elseif range_parameter2 == "deaths_off"
        p_folder = "/death-off_" 
    elseif range_parameter2 == "deaths_on"
        p_folder = "/death-on_" 
    elseif range_parameter2 == "switchings_log"
        p_folder = "/switch_" 
    end
    for j2 = 1:J2
        infer_mutation_rates(p, m, p_folder=p_folder*"$j2") 
    end
end

# Reproduce all data in the manuscript

function data_inference_manuscript()
    # Parameter regime: mutation-rate increase x switching rate
    # Estimation method: heterogeneous-response model with setting the relative division rate of response-on cells to zero (known fraction of response-on subpopulation)
    simulate_fluctuation_assays("range_mu-inc", p2="range_switching", set_seed=true)
    infer_mutation_rates("range_mu-inc", "het_zero_div", "range_switching")

    # Parameter regime: death rate of response-off x -on cells, for switching rates 0.01 and 0.05
    # Estimation method: heterogeneous-response model with setting the relative division rate of response-on cells to zero (known fraction of response-on subpopulation)
    for i in [1, 5]
        simulate_fluctuation_assays("range_death-off_switch-$i", p2="range_death-on_switch-$i", set_seed=true)
        infer_mutation_rates("range_death-off_switch-$i", "het_zero_div", "range_death-on_switch-$i")
    end

    # Parameter regime: differential fitness of response-off mutants
    # Estimation method: heterogeneous-response model with setting the relative division rate of response-on cells to zero (known fraction of response-on subpopulation)
    simulate_fluctuation_assays("range_fit-mut", set_seed=true)
    infer_mutation_rates("range_fit-mut", "het_zero_div")

    # Parameter regime: relative division rate of response-on cells
    # Estimation methods
    # (i) Heterogeneous-response model with setting the relative division rate of response-on cells to zero/true value or inferring it (known fraction of response-on subpopulation)
    # (ii) Heterogeneous-response model with setting the relative division rate of response-on cells to zero (unknown fraction of response-on subpopulation)
    # (iii) Homogeneous-response model without/with/jointly inferring the differential fitness of mutants
    simulate_fluctuation_assays("range_rel-div-on", set_seed=true)
    for m in ["het_zero_div", "het_set_div", "het_infer_div", "het_zero_div_unknown_fraction", "hom_no_fit", "hom_infer_fit", "hom_joint_fit"]
        infer_mutation_rates("range_rel-div-on", m)
    end
end

function data_supplementary_material()
    # Parameter regime: switching rate x relative division rate of response-on cells
    # Simulating response-on non-mutants stochastically to test assumption A1
    simulate_fluctuation_assays("range_switching", p2="range_rel-div-on", set_seed=true, A1=true)
end
