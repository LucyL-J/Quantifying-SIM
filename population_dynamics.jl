using Roots, Distributions, SpecialFunctions

# Population growth dynamics of non-mutants (deterministic)
# Growth of response-off subpopulation: exponential growth with population growth rate = division-death-switching rate, initial population size N0
pop_size(t, N0, pop_growth) = N0*exp(pop_growth*t)
# Growth of response-on subpopulation: influx from response-off subpopulation plus exponential growth with own net growth rate = division-death rate, initial population size N0_on 
# Initial population size is set to zero
function pop_size(t, N0_off, pop_growth, switching, N0_on, net_growth_on)
    if net_growth_on == pop_growth
        return N0_off*switching*t*exp(pop_growth*t)
    else
        return N0_off*switching/(pop_growth - net_growth_on) * (exp(pop_growth*t) - exp(net_growth_on*t)) + N0_on*exp(net_growth_on*t)
    end
end

# Time-inhomogeneous Poisson process: the time until the next event can be calculated by drawing a random variable y~U(0,1) as the input for the functions below
# Time until next event (either switching or mutation) for response-off subpopulation
interval_event(t, N0, pop_growth, event_rate, y) = log(y*pop_growth/(event_rate*pop_size(t, N0, pop_growth)) + 1) / pop_growth
# Time until next mutation event for response-on subpopulation
function interval_event(t, N0_off, pop_growth, switching, N0_on, net_growth_on, T, event_rate, y)
    if net_growth_on == 0.
        f1(tau) = event_rate * (N0_off*switching/pop_growth * (exp(pop_growth*(t + tau))/pop_growth - exp(pop_growth*t)/pop_growth) - (N0_off*switching/pop_growth - N0_on) * tau) - y
        try
            return find_zero(f1, (0., T))
        catch e
            return T
        end
    elseif net_growth_on == pop_growth
        f2(tau) = event_rate/pop_growth * ((N0_off*switching*tau+N0_on) * exp(pop_growth*(t+tau)) + N0_off*switching*(t-1/pop_growth)/pop_growth * (exp(pop_growth*tau) - 1) - N0_on*exp(pop_growth*tau)) - y
        try
            return find_zero(f2, (0., T))
        catch e
            return T
        end
    else
        f3(tau) = event_rate * (N0_off*switching/(pop_growth-net_growth_on) * (exp(pop_growth*(t + tau))/pop_growth - exp(pop_growth*t)/pop_growth) - (N0_off*switching/(pop_growth-net_growth_on) - N0_on) * (exp(net_growth_on*(t + tau))/net_growth_on - exp(net_growth_on*t)/net_growth_on)) - y
        try
            return find_zero(f3, (0., T))
        catch e
            return T
        end
    end
end
# Expected time of the first mutation
t_first_m(N0, pop_growth, mutation) = gamma(0., mutation*N0/pop_growth) / pop_growth
# Expected time of reaching a final population size Nf
t_final(N0, pop_growth, Nf) = log(Nf/N0) / pop_growth
# Time when the expected number of mutations (in the total population) is given by M; this is used to determine the duration of the growth phase of the simulated fluctuation assays
function t_expected_m(N0_off, pop_growth, mutation_off, switching, N0_on, net_growth_on, mutation_on, M)
    if net_growth_on == 0.
        f1(t) = N0_off/pop_growth*(mutation_off + switching*mutation_on/pop_growth)*(exp(pop_growth*t) - 1) + mutation_on*(N0_on - switching*N0_off/pop_growth)*t - M
        return find_zero(f1, 0.)
    elseif net_growth_on == pop_growth
        f2(t) = (mutation_off*N0_off/pop_growth + mutation_on*(N0_on - switching*N0_off/pop_growth)/pop_growth) * (exp(pop_growth*t) - 1) + mutation_on*switching*N0_off*t*exp(pop_growth*t)/pop_growth - M
        return find_zero(f2, 0.)
    else
        f3(t) = N0_off*(mutation_off + switching*mutation_on/(pop_growth-net_growth_on))*(exp(pop_growth*t) - 1)/pop_growth + mutation_on*(N0_on - switching*N0_off/(pop_growth-net_growth_on))*(exp(net_growth_on*t) - 1)/net_growth_on - M
        return find_zero(f3, 0.)
    end
end

# Branching processes (single/two-type)
# Single-type branching process
function number_cells(T, birth, death; N_max=Inf)
    if T <= 0.
        return 0
    else
        if birth == 0. && death == 0.
            return 1
        else
            m = 1
            t = rand(Exponential(1/((birth+death)*m)))
            while t < T
                if m >= N_max
                    m = Int(round(pop_size(T-t, m, birth-death)))
                    break
                end
                r = rand()
                if r <= birth/(birth+death)
                    m += 1
                else
                    m -= 1
                end
                if m == 0
                    return 0
                    break
                end
                tau = rand(Exponential(1/((birth+death)*m)))
                t += tau
            end
            return m
        end
    end
end
# Two-type branching process with switching from response-off to response-on
function number_cells(T, birth_off, death_off, switching, birth_on, death_on)
    if T <= 0.
        return [0, 0]
    else
        m_off = 1
        m_on = 0
        t = rand(Exponential(1/((birth_off+death_off+switching)*m_off)))
        while t < T
            r = rand()
            if r <= birth_off/(birth_off+death_off+switching)
                m_off += 1
            elseif r <= (birth_off+switching)/(birth_off+death_off+switching)
                m_off -= 1
                m_on += number_cells(T-t, birth_on, death_on)
            else
                m_off -= 1
            end
            if m_off == 0
                return [0, m_on]
                break
            end
            tau = rand(Exponential(1/((birth_off+death_off+switching)*m_off)))
            t += tau
        end
        return [m_off, m_on]
    end
end

# Simulating the population dynamics during the growth phase of a fluctuation assay
# Mutant count for a homogeneous population with cell death and differential fitness of mutants
function mutant_count(T, N0, division, death, fitness_m, mutation)
    t = 0.
    m = 0
    while t < T
        y = rand(Exponential(1))
        tau = interval_event(t, N0, division-death, mutation, y)
        t += tau
        m += number_cells(T-t, division*fitness_m, death)
    end
    return m
end
# Mutant count data for a homogeneous population and a given number of cultures
# Optional: cell death and differential fitness of mutants
# Returns mutant counts and final population size
function mutant_count(T, N0, division, mutation, num_cultures::Int; death=0., fitness_m=1.)
    mc = Vector{Int}(undef, num_cultures)
    for i = 1:num_cultures
        mc[i] = mutant_count(T, N0, division, death, fitness_m, mutation)
    end
    Nf = pop_size(T, N0, division-death)
    return mc, Nf
end
# Mutant count for a heterogeneous population with response-off and -on subpopulation (division, death and mutation rate each and switching from response-off to -on)
function mutant_count(T, N0_off, division_off, death_off, fitness_m_off, mutation_off, switching, N0_on, division_on, death_on, mutation_on)
    pop_growth = division_off - death_off - switching
    t = 0.
    m = 0
    while t < T
        y = rand(Exponential(1))
        tau = interval_event(t, N0_off, pop_growth, mutation_off, y)
        t += tau
        m += sum(number_cells(T-t, division_off*fitness_m_off, death_off, switching, division_on, death_on))
    end
    t = 0.
    while t < T
        y = rand(Exponential(1))
        tau = interval_event(t, N0_off, pop_growth, switching, N0_on, division_on-death_on, T, mutation_on, y)
        t += tau
        m += number_cells(T-t, division_on, death_on)
    end
    return m
end
# Mutant count data for a heterogeneous population and a given number of cultures 
# Optional: cell death of response-off and/or -on cells, differential fitness of response-off mutants and initial fraction of response-on subpopulation
# Returns mutant counts, final population size and final fraction of response-on subpopulation 
function mutant_count(T, N0, division_off, mutation_off, switching, division_on, mutation_on, num_cultures::Int; death_off=0., fitness_m_off=1., death_on=0., f0_on=switching/(division_off-death_off))
    mc = Vector{Int}(undef, num_cultures)
    for i = 1:num_cultures
        mc[i] = mutant_count(T, N0*(1-f0_on), division_off, death_off, fitness_m_off, mutation_off, switching, N0*f0_on, division_on, death_on, mutation_on)
    end
    Nf_off = pop_size(T, N0*(1-f0_on), division_off-death_off-switching)
    Nf_on = pop_size(T, N0*(1-f0_on), division_off-death_off-switching, switching, N0*f0_on, division_on-death_on)
    return mc, Nf_off+Nf_on, Nf_on/(Nf_off+Nf_on)
end



# Testing validity of approximations
# Simulating non-mutant dynamics: response-off subpopulation deterministic and response-on subpopulation stochastic
# Returns the number of response-on cells at time T
function non_mutants_on(T, N0_off, pop_growth, switching, N0_on::Int, division_on, death_on, N_max)
    if N0_on >= N_max
        return Int(round(pop_size(T, N0_off, pop_growth, switching, N0_on, division_on-death_on)))
    else
        t = 0.
        n_on = 0
        for i = 1:N0_on
            n_on += number_cells(T, division_on, death_on, N_max=N_max)
        end
        while t < T
            y = rand(Exponential(1))
            tau = interval_event(t, N0_off, pop_growth, switching, y)
            t += tau
            n_on += number_cells(T-t, division_on, death_on, N_max=N_max)
        end
        return n_on
    end
end
# Simulating mutant dynamics; non-mutants deterministic and mutants (with switching) stochastic 
# Returns the number of mutants which 
# (i) come from the response-off subpopulation and did not switch from response-off to response-on
# (ii) come from the response-off subpopulation and switched from response-off to response-on
# (iii) come from the response-on subpopulation
function mutants_off_switch_on(T, N0_off, division_off, death_off, mutation_off, switching, N0_on, division_on, death_on, mutation_on)
    pop_growth = division_off - death_off - switching
    t = 0.
    m_off_switch = [0, 0]
    while t < T
        y = rand(Exponential(1))
        tau = interval_event(t, N0_off, pop_growth, mutation_off, y)
        t += tau
        m_off_switch .+= number_cells(T-t, division_off, death_off, switching, division_on, death_on)
    end
    t = 0.
    m_on = 0
    while t < T
        y = rand(Exponential(1))
        tau = interval_event(t, N0_off, pop_growth, switching, N0_on, division_on-death_on, T, mutation_on, y)
        t += tau
        m_on += number_cells(T-t, division_on, death_on)
    end
    return [m_off_switch; m_on]
end