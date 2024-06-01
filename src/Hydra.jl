module Hydra

    using Reexport
    @reexport using DEBParamStructs
    @reexport using DoseResponse
    using ComponentArrays, StaticArrays
    using Parameters, OrdinaryDiffEq
    using Distributions, StatsBase, Random
    using DataFrames
    using PrecompileTools

    abstract type AbstractSpeciesParams <: AbstractParams end
    abstract type AbstractABM end
    abstract type AbstractAgent end

    include("Solvers.jl")

    include("Params.jl")
    export childstruct!, AbstractABM, AbstractSpeciesParams, ABM, DEBAgent, GlobalParams, GlobalBaseStatevars, SpeciesParams, DEBParamCollection, AgentParams

    include("DEBABM.jl")

    include("StateVars.jl")
    export init_substates_agent, init_substates_global, initialize_statevars, initialize_statevars!, initialize_agents!

    include("IO.jl")
    export setproperty!, isolate_pmoas!, set_equal!, relative_response

    include("ModelFunctions.jl")
    export sig, clipneg

    include("DEBODE.jl")
    export abstractsimulator, returntypes, simulator, @replicates

    include("ImpliedTraits.jl")
    include("Macros.jl")

end # module Hydra
