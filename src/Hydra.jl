module Hydra

    using ComponentArrays, StaticArrays
    using Parameters, OrdinaryDiffEq
    using Distributions, StatsBase, Random
    using DataFrames
    using PrecompileTools
    using RecipesBase


    # establishing type hierarchy
    abstract type AbstractParams end
    abstract type AbstractParamCollection end # an AbstractParamCollection contain a defined set of multiple AbstractParams instances
    abstract type AbstractSpeciesParams <: AbstractParams end
    abstract type AbstractGlobalParams <: AbstractParams end

    include("DoseResponse.jl") # collection of useful dose-response functions
    export LL2, LL2h, LL2M, LL2inv, LL2hinv, WB2, WB2, LLBP5, LLAS3, LL3, CRS6, CRS4, CRS4U, CRS6U, CRS5US, NEC2pos, NEC2neg
    
    include("Params.jl") # parameter structures
    export AbstractParams, AbstractParamCollection, childstruct!, AbstractABM, AbstractSpeciesParams, ABM, GlobalParams, GlobalBaseStatevars, SpeciesParams, DEBParamCollection, AgentParams

    include("StateVars.jl") # initialization of state variables
    export init_substates_agent, init_substates_global, initialize_statevars, initialize_statevars!, initialize_agents!

    include("IO.jl") # input/output handling
    export setproperty!, isolate_pmoas!, set_equal!, relative_response

    include("ModelFunctions.jl") # core model functions (derivatives and rules)
    export sig, clipneg

    include("HydraODE.jl") # wrappers for ODE-based simulations
    export abstractsimulator, returntypes, simulator, @replicates

    include("ImpliedTraits.jl") # calculation of traits from parameters
    include("Macros.jl") # quality of life-stuff

    include("Recipes.jl") # plotting
    export lineplot, groupedlineplot, rugplot, gridxlabel, gridylabel

    
    @compile_workload begin
        theta = DEBParamCollection()
        yhat = simulator(theta)
    end

end # module Hydra
