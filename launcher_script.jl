#import Pkg
#Pkg.instantiate()
#Pkg.activate("packages")
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
                mc, Nf, f_on = mutant_count(T, N0, division_off, mutation_off_s, switching, division_on, mutation_on, num_cultures, fitness_m_off=fitness_m_s, death_off=death_off, death_on=death_on, f0_on=f0_on)
                mutant_counts[:,"$j"] = mc
                p_final[:,"$j"] = [Nf, f_on, T]
            end
            j += 1
        end
        if S1 == true
            CSV.write("output_data/"*range_2*"n_t1-"*range*"-$(r2_parameter)_$j2.csv", n_t1)
            CSV.write("output_data/"*range_2*"n_tf-"*range*"-$(r2_parameter)_$j2.csv", n_tf)
            j2 += 1
        else
            T_p = t_expected_m(N0, division_off, mutation_off_p, 0., 0, 0., 0., expected_M)
            mc_p, Nf_p = mutant_count(T_p, N0, division_off, mutation_off_p, num_cultures, fitness_m=fitness_m_p)
            mutant_counts.p = mc_p
            p_final.p = [Nf_p, 0., T_p]
            if range_2 == "" 
                CSV.write("output_data/mutant_counts-"*range*".csv", mutant_counts)
                CSV.write("output_data/p_final-"*range*".csv", p_final)
            else
                CSV.write("output_data/"*range_2*"mutant_counts-"*range*"-$(r2_parameter)_$j2.csv", mutant_counts)
                CSV.write("output_data/"*range_2*"p_final-"*range*"-$(r2_parameter)_$j2.csv", p_final)
                j2 += 1
            end
        end
    end
end

# Mutation rate estimation and model for ranges of parameters
function infer_mutation_rates(range, mod, num_cultures, range_2=""; conf=false)
    mc_bound = 1000                                                # Discarding all mutant counts >mc_bound
    R = 100
    p = DataFrame(CSV.File("input_parameters/"*range*".csv"))
    r_parameter = p.r_parameter[1]                              
    r_start = p.r_start[1]
    r_end = p.r_end[1]
    r_i = p.r_increment[1]
    if mod == "model_selection"
        est_res_1 = DataFrame()
        est_res_2 = DataFrame()
        est_res_3 = DataFrame()
        est_res_4 = DataFrame()
        est_res_5 = DataFrame()
        selected_model = DataFrame()
    else
        est_res = DataFrame()
    end
    try
        mkdir("inferred_parameters/"*range_2)
    catch e
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
        j2 = 1  
        range_2 *= "/"
    end    
    for r2 = r2_start:r2_i:r2_end
        if range_2 == "" 
            mutant_counts = DataFrame(CSV.File("output_data/mutant_counts-"*range*".csv"))
            p_final = DataFrame(CSV.File("output_data/p_final-"*range*".csv"))
            suffix = "/"
        else
            mutant_counts = DataFrame(CSV.File("output_data/"*range_2*"mutant_counts-"*range*"-$(r2_parameter)_$j2.csv"))
            p_final = DataFrame(CSV.File("output_data/"*range_2*"p_final-"*range*"-$(r2_parameter)_$j2.csv"))
            suffix = "-$(r2_parameter)_$j2/"
            j2 += 1
        end
        mc_p = mutant_counts[:,end]
        mc_p = mc_p[mc_p .< mc_bound]
        Nf_p = p_final[1,end]
        j = 1
        for r = r_start:r_i:r_end
            try
                mkdir("inferred_parameters/"*range_2*range*suffix)
            catch e
            end
            mc_s = mutant_counts[:,j]
            mc_s = mc_s[mc_s .< mc_bound]
            Nf_s = p_final[1,j]
            f_on = p_final[2,j]
            n = 1
            for c in num_cultures
                for i = 1:R
                    while (sum(mc_p[n+(i-1)*c:n-1+i*c]) == 0) || (sum(mc_s[n+(i-1)*c:n-1+i*c]) == 0)
                        n += c
                    end
                    if mod == "het_zero-div"
                        est_res[:,"$i"] = estimu_het(mc_p[n+(i-1)*c:n-1+i*c], Nf_p, mc_s[n+(i-1)*c:n-1+i*c], Nf_s, f_on, conf=conf)
                    elseif mod == "het_set-div"
                        est_res[:,"$i"] = estimu_het(mc_p[n+(i-1)*c:n-1+i*c], Nf_p, mc_s[n+(i-1)*c:n-1+i*c], Nf_s, f_on, r, conf=conf)
                    elseif mod == "het_infer-div"
                        est_res[:,"$i"] = estimu_het(mc_p[n+(i-1)*c:n-1+i*c], Nf_p, mc_s[n+(i-1)*c:n-1+i*c], Nf_s, f_on, false, conf=conf)
                    elseif mod == "model_selection"
                        hom_1, hom_2, hom_3, het_4, het_5, selected_m = estimu_select(mc_p[n+(i-1)*c:n-1+i*c], Nf_p, mc_s[n+(i-1)*c:n-1+i*c], Nf_s, conf=conf)
                        est_res_1[:,"$i"] = hom_1
                        est_res_2[:,"$i"] = hom_2
                        est_res_3[:,"$i"] = hom_3
                        est_res_4[:,"$i"] = het_4
                        est_res_5[:,"$i"] = het_5
                        selected_model[:, "$i"] = selected_m
                    end
                end
                if mod == "model_selection"
                    #selected_model.Criterion = ["LRT (hom)", "LRT (het)", "AIC", "BIC", "CV"]
                    CSV.write("inferred_parameters/"*range_2*range*suffix*"hom_wo-fitm-number_cultures_$c-$(r_parameter)_$j.csv", est_res_1)
                    CSV.write("inferred_parameters/"*range_2*range*suffix*"hom_fitm-number_cultures_$c-$(r_parameter)_$j.csv", est_res_2)
                    CSV.write("inferred_parameters/"*range_2*range*suffix*"hom_fitm-unconstr-number_cultures_$c-$(r_parameter)_$j.csv", est_res_3)
                    CSV.write("inferred_parameters/"*range_2*range*suffix*"het_zero-div_unknown-f-number_cultures_$c-$(r_parameter)_$j.csv", est_res_4)
                    CSV.write("inferred_parameters/"*range_2*range*suffix*"het_infer-div_unknown-f-number_cultures_$c-$(r_parameter)_$j.csv", est_res_5)
                    CSV.write("inferred_parameters/"*range_2*range*suffix*"selected_model-number_cultures_$c-$(r_parameter)_$j.csv", selected_model)
                else
                    CSV.write("inferred_parameters/"*range_2*range*suffix*mod*"-number_cultures_$c-$(r_parameter)_$j.csv", est_res)
                end
                n += R*c
            end
            j += 1
        end
    end
end

# Reproduce all data in the manuscript

function data_inference_manuscript()
    # Parameter regime: mutation-rate increase x relative switching rate
    # Estimation method: heterogeneous-response model with setting the relative division rate of response-on cells to zero (known fraction of response-on subpopulation)
    simulate_fluctuation_assays("range-nu_on", "range-alpha-f0", set_seed=true)
    infer_mutation_rates("range-nu_on", "het_zero-div", [50,20,10], "range-alpha-f0", conf=true)

    # Parameter regime: death rate of response-off x -on cells, for switching rates 0.01 and 0.05
    # Estimation method: heterogeneous-response model with setting the relative division rate of response-on cells to zero (known fraction of response-on subpopulation)
    for i in [1, 5]
        simulate_fluctuation_assays("range-delta_off-alpha$i", "range-delta_on-alpha$i", set_seed=true)
        infer_mutation_rates("range-delta_off-alpha$i", "het_zero-div", [50], "range-delta_on-alpha$i")
    end

    # Parameter regime: differential fitness of response-off mutants
    # Estimation method: heterogeneous-response model with setting the relative division rate of response-on cells to zero (known fraction of response-on subpopulation)
    simulate_fluctuation_assays("range-rho", set_seed=true)
    infer_mutation_rates("range-rho", "het_zero-div", [50])

    # Parameter regime: relative division rate of response-on cells
    # Estimation methods
    # (i) Heterogeneous-response model with setting the relative division rate of response-on cells to zero/true value or inferring it (known fraction of response-on subpopulation)
    # (ii) Heterogeneous-response model with setting the relative division rate of response-on cells to zero or inferring it (unknown fraction of response-on subpopulation)
    # (iii) Homogeneous-response model without/with/jointly inferring the differential fitness of mutants
    for i in [10, 100]
        simulate_fluctuation_assays("range-gamma_on-increase$i", set_seed=true)
    end
    infer_mutation_rates("range-gamma_on-increase10", "model_selection", [50], conf=true)
    infer_mutation_rates("range-gamma_on-increase100", "model_selection", [50,20,10], conf=true)
    for mod in ["het_zero-div", "het_set-div", "het_infer-div"]
        infer_mutation_rates("range-gamma_on-increase100", mod, [50], conf=true)
    end

    # Parameter regime: mutation rate of response-off cells x differential fitness of response-off mutants (homogeneous response)
    # (i) Heterogeneous-response model with setting the relative division rate of response-on cells to zero or inferring it (unknown fraction of response-on subpopulation)
    # (ii) Homogeneous-response model without/with/jointly inferring the differential fitness of mutants
    for i in ["", "_unconstr"]
        simulate_fluctuation_assays("range-nu_off_s", "range-rho"*i, set_seed=true)
        infer_mutation_rates("range-nu_off_s", "model_selection", [50], "range-rho"*i, conf=true)
    end
end

function data_supplementary_material()
    # Parameter regime: switching rate x relative division rate of response-on cells
    # Simulating response-on non-mutants stochastically to test assumption S1
    for f in ["f0", "fstat"]
        simulate_fluctuation_assays("range-gamma_on-increase100", "range-alpha-"*f, set_seed=true, S1=true)
    end
end