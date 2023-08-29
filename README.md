# Quantifying-SIM
Documentation of the simulations and computational analysis of the manuscript \
*Quantifying mutation rates under heterogeneous stress responses* 

The simulations and the inference algorithms are written in the programming language [julia](https://julialang.org). \
To replicate the data presented in the manuscript and the additional data in the supplementary material, clone this repository, install [julia verison 1.6](https://julialang.org/downloads/#long_term_support_release), start it and execute
```
cd("...")
include("launcher_script.jl")
data_inference_manuscript()
data_supplementary_material()
```
with the dots replaced by the path to the repository.

To analyise mutant count data from a pair of fluctuation assays under permissive/stressful conditions, use the script `inference.jl`. To set up, start julia and execute 
```
cd("...")
import Pkg
Pkg.activate("packages")
include("inference.jl")
```
with the dots replaced again by the path to the repository.

The function 
```
estimate_mu(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s; fitm_p=false, fitm_s=false, joint=false)
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
* `joint`: To constrain the mutant fitness to be equal under permissive and stressful conditions, set `joint=true`

As output, the function prints whether/which model is selected (difference in AIC > 2) and the respective inferred parameters.

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

Optionally, the input parameter `fit_m` can be set to a value measured in a separate experiment, as a joint inference parameter via `rel_div_on="joint"` or as two inference parameters (permissive/stressful condition) via `rel_div_on="infer"`. \
The output parameters are 
* Mutation rate permissive condition 
* Mutation rate stressful condition
* Optional: Inferred relative fitness of mutants compared to non-mutants (under permissive/stressful condition)
* AIC

From the inferred mutation rates, the increase in population mean mutation rate due to stress can be calculated.\
If the inference fails, AIC=`Inf` is returned.