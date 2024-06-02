
using Pkg; Pkg.activate("test")
using BenchmarkTools
using OrdinaryDiffEq
using Test
using Revise
using Hydra

@benchmark out = simulator(ODEParamCollection()) 
@test true # this just has to run without an error