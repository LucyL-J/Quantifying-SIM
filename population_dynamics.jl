using Roots

# Population growth dynamics of non-mutants (deterministic)
# Growth of response-off subpopulation: exponential growth with population growth rate = division-death-switching rate, initial population size N0
pop_size(t, N0, pop_growth) = N0*exp(pop_growth*t)
# Growth of response-on subpopulation: influx from response-off subpopulation plus exponential growth with own net growth rate = division-death rate, initial population size N0_on 
# Initial population size is set to zero
function pop_size(t, N0, pop_growth, switching, net_growth_on, N0_on)
    if net_growth_on == pop_growth
        return N0*switching*t*exp(pop_growth*t)
    else
        return N0*switching/(pop_growth - net_growth_on) * (exp(pop_growth*t) - exp(net_growth_on*t)) + N0_on*exp(net_growth_on*t)
    end
end

# Time-inhomogeneous Poisson process: the time until the next event can be calculated by drawing a random variable y~U(0,1) as the input for the functions below
# Time until next event (either switching or mutation) for response-off subpopulation
interval_event(t, N0, pop_growth, event_rate, y) = log(y*pop_growth/(event_rate*pop_size(t, N0, pop_growth)) + 1) / pop_growth
# Time until next mutation event for response-on subpopulation
function interval_event(t, N0, pop_growth, switching, net_growth_on, N0_on, T, event_rate, y)
    if net_growth_on == 0.
        return interval_event(t, switching*N0/pop_growth, pop_growth, event_rate, y)
    elseif net_growth_on == pop_growth
        f1(tau) = event_rate/pop_growth * ((N0*switching*tau+N0_on) * exp(pop_growth*(t+tau)) + N0*switching*(t-1/pop_growth)/pop_growth * (exp(pop_growth*tau) - 1) - N0_on*exp(pop_growth*tau)) - y
        try
            return find_zero(f1, (0., T))
        catch e
            return T
        end
    else
        f2(tau) = event_rate * (N0*switching/(pop_growth-net_growth_on) * (exp(pop_growth*(t + tau))/pop_growth - exp(pop_growth*t)/pop_growth) - (N0*switching/(pop_growth-net_growth_on) - N0_on) * (exp(net_growth_on*(t + tau))/net_growth_on - exp(net_growth_on*t)/net_growth_on)) - y
        try
            return find_zero(f2, (0., T))
        catch e
            return T
        end
    end
end