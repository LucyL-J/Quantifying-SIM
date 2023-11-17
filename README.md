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

### Estimation using the heterogeneous-response model

We developed a population dynamic model that considers cell-to-cell heterogeneity in stress responses in the form of a cell subpopulation with highly expressed stress response (compared to the rest of the population with a low expression level). The following function estimates mutation rates under this _heterogeneous-response_ model from fluctuation assay data
```
estimate_mu_het(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s, f_on; rel_div_on=0.)
```

The input parameters are
* `mc_p`: Mutant counts for permissive condition (a vector with integers)
* `Nf_p`: Avergae final population size for permissive condition
* `mc_s`: Mutant counts for stressful condition (a vector with integers)
* `Nf_s`: Average final population size for stressful condition 

Optional
* `f_on`: Fraction of the response-on subpopulation (when known from a separate experimental measurement). Inferred when not given as an iput parameter.
* `rel_div_on`: Relative division rate of response-on cells can be set to a value measured in a separate experiment or set as an additional inference parameter via `rel_div_on="infer"`. 

The output parameters are 
* Mutation rate response-off
* Mutation rate response-on
* Optional: Inferred fraction of the response-on subpopulation
* Relative fitness of response-on compared to response-off cells
* AIC

From the inferred mutation rates, the relative mutation-rate increase associated with the induction of the stress response can be calculated.\
If the inference fails, AIC=`Inf` is returned.

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

#### Case 1: The fraction of the response-on subpopulation is unknown

Evaluating the estimation function (by default, the relative division rate of response-on cells is set to zero)
```
res = estimate_mu_het(mc_p, 10^8, mc_s, 1.6*10^7)
```
yields the following output
```
5-element Vector{Float64}:
   1.0063765011105267e-8
   1.7818401948380247e-6
   0.023222420458462154
   0.0
 406.0980895989688
```
with an estimated mutation-rate increase of $\approx 177$-fold and fraction of the response-on subpopulation of $\approx 2.3 \%$. \
Instead of setting the relative division rate of response-on cells to zero, we can set it as an inference parameter
```
res = estimate_mu_het(mc_p, 10^8, mc_s, 1.6*10^7, rel_div_on="infer")
```
yielding
```
5-element Vector{Float64}:
   9.93029230085328e-9
   1.5722757467417157e-6
   0.02524485797435025
   0.0809156313118107
 407.72583389687406
```
i.e. an estimated mutation-rate increase of $\approx 158$-fold, fraction of response-on subpopulation of $\approx 2.5 \%$ and relative division rate of response-on cells of $\approx 0.08$. However, this model version fits the data less well, as it's AIC is higher.

#### Case 2: An estimate of the fraction of the response-on subpopulation is available from a separate experiment

Let's consider the case that the fraction of cells with elevated stress response has been estimated to be around $\approx 5\%$. Then, evalutating the estimation function
```
res = estimate_mu_het(mc_p, 10^8, mc_s, 1.6*10^7, 0.05)
```
yields
```
4-element Vector{Float64}:
   1.0063709988215868e-8
   8.048868057992385e-7
   0.0
 404.0980895994152
```
with an estimated mutation-rate increase of $\approx 80$-fold. \
Again, the relative division rate of response-on cells can be inferred via setting `rel_div_on="infer`.

### Estimation using the homogeneous-response model

We also implement the standard method using a _homogeneous-response_ model with optional differential mutant fitness. A novel estimation option is to constrain the mutant fitness to be equal under permissive and stressful conditions. 

The function 
```
estimate_mu_hom(mc_p::Vector{Int}, Nf_p, mc_s::Vector{Int}, Nf_s; fit_m=1.)
```
estimates population-wide mutation rates from fluctuation assay data. \
The input parameters are
* `mc_p`: Mutant counts for permissive condition (a vector with integers)
* `Nf_p`: Average final population size for permissive condition
* `mc_s`: Mutant counts for stressful condition (a vector with integers)
* `Nf_s`: Average final population size for stressful condition

Optionally, the mutant fitness `fit_m` can be set to a value measured in a separate experiment, set as a joint inference parameter via `fit_m="joint"` or as two separate inference parameters (permissive/stressful condition) via `fit_m="infer"`. 

The output parameters are 
* Mutation rate permissive condition 
* Mutant fitness permissive condition (relative to non-mutants)
* Mutation rate stressful condition
* Mutant fitness stressful condition (relative to non-mutants)
* AIC

From the inferred mutation rates, the increase in population mean mutation rate due to stress can be calculated.\
If the inference fails, AIC=`Inf` is returned.

**Example (continued)** 

Evaluating the estimation function (by default, the mutant fitness is set to 1)
```
res = estimate_mu_hom(mc_p, 10^8, mc_s, 1.6*10^7)
```
yields the following output
```
5-element Vector{Float64}:
   1.036812895645713e-8
   1.0
   3.973310967312984e-8
   1.0
 418.87316975292185
```
and an estimated increase in population mean mutation rate of $\approx 3.8$-fold.

We can also infer the mutant fitness. Either as two separate inference parameters under permissive/stressful condtitions via 
```
res = estimate_mu_hom(mc_p, 10^8, mc_s, 1.6*10^7, fit_m="infer")
```
which yields
```
5-element Vector{Float64}:
   9.892353672581243e-9
   1.2525855391697718
   4.9646745015561864e-8
   0.13774260645934688
 403.32920230179843
```
i.e. an estimated increase in population mean mutation rate of $\approx 5$-fold and mutant fitness under permissive/stressful consitions of $\approx 1.25$ and $\approx 0.14$, respectively. \
Alternatively, we can constrain the mutant fitness to be equal under permissive/stressful conditions via
```
res = estimate_mu_hom(mc_p, 10^8, mc_s, 1.6*10^7, fit_m="joint")
```
yielding 
```
5-element Vector{Float64}:
   1.0628312661875e-8
   0.8882119940565083
   4.019999797782713e-8
   0.8882119940565083
 420.2910160915429
```
with an estimated increase in population mean mutation rate of $\approx 3.8$-fold and a mutant fitness of $\approx 0.89$.