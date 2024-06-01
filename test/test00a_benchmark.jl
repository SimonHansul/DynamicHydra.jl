
using Pkg; Pkg.activate("test")
using BenchmarkTools
using OrdinaryDiffEq
using Test
using Revise
using Hydra

@benchmark out = simulator(DEBParamCollection()) 
@test true # this just has to run without an error