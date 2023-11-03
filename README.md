# Quantifying-SIM
Documentation of the simulations and computational analysis of the manuscript \
*Estimating mutation rates under heterogeneous stress responses* \
Lucy Lansch-Justen, Meriem El Karoui, Helen K. Alexander \
[bioRxiv 2023.09.05.555499](https://doi.org/10.1101/2023.09.05.555499).

The simulations and the inference algorithms are written in the programming language [julia verison 1.6](https://julialang.org/downloads/#long_term_support_release). 

## Reproducing the data in the manuscript
To replicate the data presented in the manuscript, clone this repository, start julia and execute
```
cd("...")
include("launcher_script.jl")
data_inference_manuscript()
data_supplementary_material()
```
with the dots replaced by the path to the repository. \
The raw figures from the manuscript can be generated using the notebook `plots.ipynb`.

## Analysing experimental mutant count data

To analyise mutant count data from a pair of fluctuation assays under permissive/stressful conditions, use the script `inference.jl`. To set up, start julia and execute 
```
cd("...")
import Pkg
Pkg.activate("packages")
Pkg.instantiate()
include("inference.jl")
```
with the dots replaced again by the path to the repository.

### Case 1: It is not known whether the stress response is heterogeneous and/or there is no estimate of the fraction of cells with elevated expression level of the response available.

The function 
```
estimate_mu(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s; fitm_p=1., fitm_s=1., infer_fitm=true, fit_on=0., infer_fit_on=true)
``` 
compares between the heterogeneous- and the homogeneous-response model using the AIC and estimates the mutation rates under the selected model. \
The input parameters are
* `mc_p`: Mutant counts for permissive condition (has to be a vector with integers)
* `Nf_p`: Final population size for permissive condition
* `mc_s`: Mutant counts for stressful condition (has to be a vector with integers)
* `Nf_s`: Final population size for stressful condition

Optional input parameters are
* `fitm_p`: Mutant fitness under permissive condition (relative to non-mutants)
* `fitm_s`: Mutant fitness under stressful condition (relative to non-mutants)
* `infer_fitm`: To constrain the mutant fitness to be equal under permissive and stressful conditions, set `infer_fitm="joint"`. To not consider mutant fitness as inference parameters, set `infer_fitm=false`
* `fit_on`: Relative fitness of response-on cells (compared to response-off cells)
* `infer_fit_on`: To not consider relative fitness of response-on cells as inference parameter, set `infer_fit_on=false`

As output, the function prints whether/which model is selected (difference in AIC > 2) and the respective inferred parameters.

**Example** 

Under permissive conditions, the mutant counts 
```
mc_p = [0, 9, 3, 11, 0, 0, 1, 0, 3, 5, 1, 1, 2, 0, 0, 1, 89, 0, 1, 1, 9, 0, 0, 0, 1, 1, 0, 0, 31, 3, 0, 10, 13, 0, 3, 3, 165, 2, 0, 1, 0, 0, 0, 2, 55, 11, 3, 5, 1, 34]
```
and a final population size of $10^8$ were observed. \
Under stressful conditions, the mutant counts 
```
mc_s = [2, 4, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 2, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 2, 1, 1, 1, 1, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 2, 2, 4, 1, 0, 4, 0, 1, 0, 2, 0]
```
and a final population size of $1.6\cdot 10^7$ were observed. 

Evaluating the estimation function, constraining the mutant fitness to be equal under stressful as under permissive conditions 
```
estimate_mu(mc_p, 10^8, mc_s, 1.6*10^7, infer_fitm="joint")
```
we get the following output:
```
Heterogeneous-response model is selected
Relative fitness of response-on cells set to input value/to 0
Mutation rate response-off cells = 1.0063765009509016e-8
Mutation rate response-on cells = 1.7818401953824812e-6
Fraction response-on subpopulation = 0.02322242044793327
Relative mutation-rate increase = 177.05502798394656
Increase in population mean mutation rate = 5.088423881815865
```
Without constraining the mutant fitness, we get
```
estimate_mu(mc_p, 10^8, mc_s, 1.6*10^7)
```
```
Homogeneous-response model is selected
Mutant fitnesses inferred
Mutation rate permissive condition = 9.892353672581243e-9
Mutation rate stressful condition = 4.9646745015561864e-8
Mutant fitness permissive condition = 1.2525855391697718
Mutant fitness stressful condition = 0.13774260645934688
Increase in population mean mutation rate = 5.018698952623212
```

### Case 2: Heterogeneous-response model: an estimate of the fraction of cells with elevated expression level of the response is available.

The function 
```
estimate_mu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, f_on; rel_div_on=0.)
```
estimates the mutation rates under the heterogeneous-response model in the case that the fraction of the response-on subpopulation (under stress) is known from a separate experimental measurement. \
The input parameters are
* `mc_p`: Mutant counts for permissive condition (has to be a vector with integers)
* `Nf_p`: Final population size for permissive condition
* `mc_s`: Mutant counts for stressful condition (has to be a vector with integers)
* `Nf_s`: Final population size for stressful condition
* `f_on`: Fraction of the response-on subpopulation 

Optionally, the input parameter `rel_div_on` can be set to a value measured in a separate experiment or as an additional inference parameter via `rel_div_on="infer"`. \
The output parameters are 
* Mutation rate response-off
* Mutation rate response-on
* Optional: Inferred relative fitness of response-on compared to response-off cells
* AIC

From the inferred mutation rates, the relative mutation-rate increase associated with the induction of the stress response can be calculated.\
If the inference fails, AIC=`Inf` is returned.

**Example (continued)** 

The fraction of cell with elevated stress response is estimated to be around $5\%$. We can estimate the mutation rates via

```
mu_off, mu_on, div_on, AIC = estimate_mu_het(mc_p, 10^8, mc_s, 1.6*10^7, 0.05)
```
and calculate the mutation-rate increase
```
mu_on/mu_off
79.97913363234322
```
We can also infer the relative division rate of cells with evelated stress response via 
```
mu_off, mu_on, div_on, AIC = estimate_mu_het(mc_p, 10^8, mc_s, 1.6*10^7, 0.05, rel_div_on="infer")
```
and calculate the mutation-rate increase for the estimated relative fitness
```
mu_on/mu_off, div_on
(78.07974840080496, 0.07977206759158645)
```

### Case 3: Homogeneous-response model

The function 
```
estimate_mu_hom(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s; fit_m=1.)
```
estimates the mutation rates under the homogeneous-response model. \
The input parameters are
* `mc_p`: Mutant counts for permissive condition (has to be a vector with integers)
* `Nf_p`: Final population size for permissive condition
* `mc_s`: Mutant counts for stressful condition (has to be a vector with integers)
* `Nf_s`: Final population size for stressful condition

Optionally, the input parameter `fit_m` can be set to a value measured in a separate experiment, as a joint inference parameter via `fit_m="joint"` or as two inference parameters (permissive/stressful condition) via `fit_m="infer"`. \
The output parameters are 
* Mutation rate permissive condition 
* Mutation rate stressful condition
* Optional: Inferred relative fitness of mutants compared to non-mutants (under permissive/stressful condition)
* AIC

From the inferred mutation rates, the increase in population mean mutation rate due to stress can be calculated.\
If the inference fails, AIC=`Inf` is returned.

**Example (continued)** 

We can also use the standard model, i.e. assuming a homogeneous stress response. By default, the mutant fitness is set to 1
```
mu_p, rho_p, mu_s, rho_s, AIC = estimate_mu_hom(mc_p, 10^8, mc_s, 1.6*10^7)
```
and we can calculate the mutation-rate increase
```
mu_s/mu_p
3.832235290114035
```
It also possible to infer the mutant fitness. Either as two separate inference parameters under permissive/stressful condtitions via 
```
mu_p, rho_p, mu_s, rho_s, AIC = estimate_mu_hom(mc_p, 10^8, mc_s, 1.6*10^7, fit_m="infer")
```
and calculate the mutation-rate increase for the estimated mutant fitnesses 
```
mu_s/mu_p, rho_p, rho_s
(5.018698952623212, 1.2525855391697718, 0.13774260645934688)
```

Or we constrain the mutant fitness to be equal under stressful as under permissive conditions
```
mu_p, rho, mu_s, rho, AIC = estimate_mu_hom(mc_p, 10^8, mc_s, 1.6*10^7, fit_m="joint")
```
and calculate the mutation-rate increase for the estimated (joint) mutant fitness
```
mu_s/mu_p, rho
(3.7823499605943303, 0.8882119940565083)
```