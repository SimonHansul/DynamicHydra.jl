# Hydra.jl

Hydra is a package for flexible and modular dynamic modelling of chemical effects. <br>
The core of the package is an implementation of a reserveless Dynamic Energy Budget (DEB) model, 
which can be used for organism-level simulations via DEBODE and population simulations via DEBABM. <br>
The initial commit provides the raw source code. Examples, tests, and documentation will follow.

## Getting started

Hydra provides a set of default parameters which can be used to simulate a reserveless Dynamic Energy Budget model. 
The parameters roughly describe the life-history of *Daphnia magna* with model currency $\mu g C$, but have no further importance except that they serve as a starting point. <br>

The following is the minimal code to run a model:
```
using Hydra

params = DEBParamCollection()
out = simulator(params)
```

The output is a DataFrame containing the model trajectories, with major state variables time `t`, structural mass `S`, maturity `H` and reproduction buffer `R`.


## Relation to similar packages

In general, the design philosophy behind `Hydra` is to facilitate modular modelling approaches. 
It should be possible to combine self-sufficient models into new models without touching (or copy-pasting) the original model code. <br>

[DynamicEnergyBudgets.jl](https://github.com/BiophysicalEcology/DynamicEnergyBudgets.jl) provides facilities for DEB modelling in Julia. <br>