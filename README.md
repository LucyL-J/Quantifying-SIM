# Quantifying-SIM
Documentation of the simulations and computational analysis of the manuscript \
*Quantifying stress-induced mutagenesis using fluctuation assays* 

The simulations and the inference algorithms are written in the programming language [julia](https://julialang.org). \
To replicate the data presented in the manuscript and the supplementary material, clone this repository, install [julia verison 1.6](https://julialang.org/downloads/#long_term_support_release), start it and execute the following
```
cd("...")
include("launcher_script.jl")
data_inference_manuscript()
data_supplementary_material()
```
with the dots replaced by the path to the repository.
