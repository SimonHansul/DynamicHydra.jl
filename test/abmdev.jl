@info("Loading packages")
using Pkg; Pkg.activate("test")
using Plots, StatsPlots, Plots.Measures
using StatsBase
using Glob
using DataFrames, DataFramesMeta
using ProgressMeter
default(leg = false, lw = 1.5)

using ComponentArrays, StaticArrays
using Test
using Parameters
using Distributions
using OrdinaryDiffEq

