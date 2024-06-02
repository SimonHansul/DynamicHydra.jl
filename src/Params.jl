
#=
## Top-level abstract types

=#

abstract type AbstractParams end
abstract type AbstractParamCollection end

copy(theta::AbstractParams) = typeof(theta)([getproperty(theta, x) for x in fieldnames(typeof(theta))]...)
copy(theta::AbstractParamCollection) = typeof(theta)([getproperty(theta, x) for x in fieldnames(typeof(theta))]...)

#=
## Constructing child structures
=#

"""
Define a childstruct
"""
function childstruct(
    parentstruct, addnames, addtypes, adddefaults
    )

    childstructname = Symbol("ChildParamType_"*randstring('a':'z', 10))

    defaultparentstruct = parentstruct() # get an instance 
    parentabstracttype = supertype(parentstruct) # we assume that the parentstruct has an  abstract supertype
    fields = vcat(fieldnames(parentstruct)..., addnames...) # get the full list of fields
    types = vcat(fieldtypes(parentstruct)..., addtypes) # get the full list of types

    defaults = vcat(
        [getproperty(defaultparentstruct, field) for field in fieldnames(parentstruct)],
        adddefaults
    )

    fielddefs = [:($(field)::$(type)=$(default)) for (field,type,default) in zip(fields,types,defaults)]

    quote 
        @with_kw mutable struct $childstructname <: $parentabstracttype
            $(fielddefs...)
        end 
    end |> eval

    return eval(quote $childstructname end)

    #@info "Defined new structure $childstructname <: $parentabstracttype from $parentstruct"
    #eval(quote $childstructname() end)
end

#=
## General structures
=#

"""
`GlobalODEParams` contain the global parameters (simulated timespan `t_max`, nutrient influx rate `Xdot_in`, etc.)
"""
@with_kw mutable struct GlobalODEParams <: AbstractParams
    N0::Int64 = 1 #  initial number of individuals [#]
    t_max::Float64 = 21. # maximum simulation time [t]
    Xdot_in::Float64 = 1200. # resource influx rate [m/t]
    k_V::Float64 = 0.1 # chemostatic dilution rate [t^-1]
    V_patch::Float64 = 0.05 # volume of a patch (or the entire similated environment) [V]
    C_W::Vector{Float64} = [0.] # external chemical concentrations [m/t], [n/t], ...
    saveat::Float64 = 1. # when to save output [t]
    odefuncs::Vector{Function} = Function[C_Wdot_const!, X_pdot_chemstat!] # ODE-based global step functions
end

"""
`GlobalODEParams` contain the global parameters (simulated timespan `t_max`, nutrient influx rate `Xdot_in`, etc.)
"""
@with_kw mutable struct GlobalABMParams <: AbstractParams
    N0::Int64 = 1 #  initial number of individuals [#]
    t_max::Float64 = 21. # maximum simulation time [t]
    Xdot_in::Float64 = 1200. # resource influx rate [m/t]
    k_V::Float64 = 0.1 # chemostatic dilution rate [t^-1]
    V_patch::Float64 = 0.05 # volume of a patch (or the entire similated environment) [V]
    C_W::Vector{Float64} = [0.] # external chemical concentrations [m/t], [n/t], ...
    AgentType::DataType = DEBAgent # the type of agent simulated
    recordagentvars::Bool = true # record agent-level output?
    saveat::Float64 = 1. # when to save output [t]
    odefuncs::Vector{Function} = Function[C_Wdot_const!, X_pdot_chemstat!] # ODE-based global step functions
    rulefuncs::Vector{Function} = Function[N_tot!] # rule-based global step functions
end



"""
`SpeciesParams` contain population means of DEB and TKTD parameters. Default values are for Daphnia magna and Azoxystrobin in μg C. <br>
Hydra.jl uses a hierarchical modelling approach where the `SpeciesParams` are  parameters which are common across all agents of a species, 
and `AgentParams` contain parameters which are specific for a species. <br>
Variability is given by the zoom factor `Z::Distribution`, which is always applied to the surface-area specific ingestion rates 
and can optionally propagate to parameters indicated in `propagate_zoom::NamedTuple`. <br>
`Z` is `Dirac(1)` by default, i.e. there is no agent variability in the default parameters. <br>
"""
@with_kw mutable struct SpeciesParams <: AbstractSpeciesParams
    Z::Distribution = Dirac(1.) # agent variability is accounted for in the zoom factor. this can be set to a Dirac distribution if a zoom factor should be applied without introducing agent variability.
    Z_male::Float64 = 1. # zoom factor for males
    sex_ratio::Float64 = 1. # initial sex ratio (females vs males)
    propagate_zoom::@NamedTuple{X_emb_int::Bool, H_p::Bool, K_X::Bool} = (X_emb_int = true, H_p = true, K_X = true) # Parameters to which Z will be propagated. Z is *always* applied to `Idot_max_rel` (with appropriate scaling).
    X_emb_int::Float64 = 19.42 # initial vitellus [m]
    K_X::Float64 = 1. # half-saturation constant for food uptake [m/V]
    Idot_max_rel::Float64 = 22.9 # maximum size-specific ingestion rate [m m^-2/3 t^-1]
    Idot_max_rel_emb::Float64 = 22.9 # size-specific embryonic ingestion rate [m m^-2/3 t^-1]
    kappa::Float64 = 0.539 # somatic allocation fraction [-]
    eta_IA::Float64 = 0.33 # assimilation efficiency [-]
    eta_AS::Float64 = 0.8 # growth efficiency [-]
    eta_SA::Float64 = 0.8 # shrinking efficiency [-]
    eta_AR::Float64 = 0.95 # reproduction efficiency [-]
    k_M::Float64 = 0.59 # somatic maintenance rate constant [t^-1]
    k_J::Float64 = 0.504 # maturity maintenance rate constant [t^-1]
    H_p::Float64 = 100. # maturity at puberty [m]
    e_S::Float64 = 0.5 # sensitivity parameter for starvation mortality (median effective S_0) [m]
    b_S::Float64 = 5. # slope parameter for starvation mortality [-]
    
    k_D_G::Vector{Float64} = [0.]  # Dominant rate constants   [t^-1] | PMoA growth efficiency
    k_D_M::Vector{Float64} = [0.]  # Dominant rate constants   [t^-1] | PMoA maintenance costs
    k_D_A::Vector{Float64} = [0.]  # Dominant rate constants   [t^-1] | PMoA assimilation efficiency
    k_D_R::Vector{Float64} = [0.38] # Dominant rate constants  [t^-1] | PMoA reproduction efficiency
    k_D_h::Vector{Float64} = [0.]   # Dominant rate constants  [t^-1] | PMoA hazard rate
    
    drc_functs_G::Vector{Function} = [LL2] # Dose-response functions | PMoA growth efficiency
    drc_functs_M::Vector{Function} = [LL2M] # Dose-response functions | PMoA maintenance costs
    drc_functs_A::Vector{Function} = [LL2] # Dose-response functions | PMoA assimilation efficiency
    drc_functs_R::Vector{Function} = [LL2] # Dose-response functions | PMoA reproduction efficiency
    drc_functs_h::Vector{Function} = [LL2h] # Dose-response functions | PMoA hazard rate

    e_G::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA growth efficiency
    e_M::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA maintenance costs
    e_A::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA assimilation efficiency
    e_R::Union{Nothing,Vector{Float64}} = [167.] # sensitivity parameters | PMoA reproduction efficiency
    e_h::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA hazard rate

    b_G::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA growth efficiency
    b_M::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA maintenance costs
    b_A::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA assimilation efficiency
    b_R::Union{Nothing,Vector{Float64}} = [0.93] # slope parameters | PMoA reproduction efficiency
    b_h::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA reproduction efficiency

    c_G::Union{Nothing,Vector{Float64}} = nothing # placeholders for additional TD parameters. only relevant if custom drc_functs with more than two parameters are defined
    c_M::Union{Nothing,Vector{Float64}} = nothing
    c_A::Union{Nothing,Vector{Float64}} = nothing
    c_R::Union{Nothing,Vector{Float64}} = nothing
    c_h::Union{Nothing,Vector{Float64}} = nothing

    d_G::Union{Nothing,Vector{Function}} = nothing
    d_M::Union{Nothing,Vector{Function}} = nothing
    d_A::Union{Nothing,Vector{Function}} = nothing
    d_R::Union{Nothing,Vector{Function}} = nothing
    d_h::Union{Nothing,Vector{Function}} = nothing

    odefuncs::Vector{Function} = Function[ # ode-based model functions
            y_z!, 
            h_S!,
            Idot!, 
            Adot!, 
            Mdot!, 
            Jdot!, 
            Sdot!, 
            Hdot!,
            H_bdot!,
            Rdot!,
            Ddot!,
            age!
        ]
    rulefuncs::Vector{Function} = Function[ # rule-based model functions
        reproduce_opportunistic!,
        die!
    ]

    adata::Vector{Symbol} = [:S, :R, :H] # statevars to record
end

"""
    AgentParams(spc::AbstractParams)

AgentParams are subject to agent variability. 
This is in contrast to SpeciesParams, which define parameters on the species-level, i.e. the population means.
"""
@with_kw mutable struct AgentParams <: AbstractParams
    Z::Float64
    Idot_max_rel::Float64
    Idot_max_rel_emb::Float64
    X_emb_int::Float64
    H_p::Float64
    K_X::Float64
    
    """
    Initialize AgentParams from SpeciesParams `spc`.
    """
    function AgentParams(spc::AbstractParams)
        agn = new()
        agent_variability!(agn, spc)
        return agn
    end
end


"""
    agent_variability!(p::Ref{AbstractParams})
Induce agent variability in spc parameters via zoom factor `Z`. 
`Z` is sampled from the corresponding distribution given in `p` and assumed to represent a ratio between maximum structurel *masses* (not lengths), 
so that the surface area-specific ingestion rate `Idot_max_rel` scales with `Z^(1/3)` and parameters which represent masses or energy pools scales with `Z`.
"""
function agent_variability!(agn::AGN, spc::SPC) where {AGN <: AbstractParams, SPC <: AbstractParams}
    agn.Z = rand(spc.Z) # sample zoom factor Z for agent from distribution
    agn.Idot_max_rel = spc.Idot_max_rel * agn.Z^(1/3) # Z is always applied to Idot_max_rel
    agn.Idot_max_rel_emb = spc.Idot_max_rel_emb * agn.Z^(1/3) #, including the value for embryos

    for param in fieldnames(typeof(spc.propagate_zoom)) # iterate over other parameters which may be affected by Z
        if getproperty(spc.propagate_zoom, param) # check whether propagation of Z should occur for this parameter
            setproperty!(agn, param, getproperty(spc, param) * agn.Z) # assign the agent value by adjusting the hyperparameter
        else # if Z should not be propagated to this parameter, 
            setproperty!(agn, param, getproperty(spc, param)) # set the agent-specific value equal to the population mean
        end
    end 
end

"""
A `ODEParamCollection` contains global parameters `glb` and spc parameters `spc`. <br>
Initialize the default parameter collection with `ODEParamCollection()`.
"""
@with_kw mutable struct ODEParamCollection <: AbstractParamCollection
    glb::AbstractParams = GlobalODEParams()
    spc::AbstractParams = SpeciesParams()
    agn::Union{Nothing,AbstractParams} = nothing
end


@with_kw mutable struct ABMParamCollection <: AbstractParamCollection
    glb::AbstractParams = GlobalABMParams()
    spc::AbstractParams = SpeciesParams()
    agn::Union{Nothing,AbstractParams} = nothing
end