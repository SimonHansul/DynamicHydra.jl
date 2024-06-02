module Hydra

    using ComponentArrays, StaticArrays
    using Parameters, OrdinaryDiffEq
    using Distributions, StatsBase, Random
    using DataFrames
    using PrecompileTools

    # establishing type hierarchy
    abstract type AbstractParams end
    abstract type AbstractParamCollection end # an AbstractParamCollection contain a defined set of multiple AbstractParams instances
    abstract type AbstractSpeciesParams <: AbstractParams end
    abstract type AbstractGlobalODEParams <: AbstractParams end
    abstract type AbstractABM end
    abstract type AbstractAgent end

    include("DoseResponse.jl") # collection of useful dose-response functions
    export LL2, LL2h, LL2M, LL2inv, LL2hinv, WB2, WB2, LLBP5, LLAS3, LL3, CRS6, CRS4, CRS4U, CRS6U, CRS5US, NEC2pos, NEC2neg
    
    include("Params.jl")
    export AbstractParams, AbstractParamCollection, childstruct!, AbstractABM, AbstractSpeciesParams, ABM, DEBAgent, GlobalODEParams, GlobalBaseStatevars, SpeciesParams, ODEParamCollection, AgentParams
    
    include("Solvers.jl") # solvers used within DEBABM
    include("DEBABM.jl") # agent-based simulations

    include("StateVars.jl") # initializeation of state variables
    export init_substates_agent, init_substates_global, initialize_statevars, initialize_statevars!, initialize_agents!

    include("IO.jl") # input/output handling
    export setproperty!, isolate_pmoas!, set_equal!, relative_response

    include("ModelFunctions.jl") # core model functions (derivatives and rules)
    export sig, clipneg

    include("DEBODE.jl") # ODE-based simulations
    export abstractsimulator, returntypes, simulator, @replicates

    include("ImpliedTraits.jl") # calculation of traits from parameters
    include("Macros.jl") # quality of life-stuff

end # module Hydra
