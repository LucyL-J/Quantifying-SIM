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