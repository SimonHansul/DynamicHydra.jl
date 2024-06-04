
# Initialize.jl
# functions to initial states, agents, models, etc. 
# note that many initializations also occur through constructors defined in Structures.jl. here we deal mainly with ComponentArrays

const X_EMB_INT_REL::Float64 = 0.001 # the (assumed) initial amount of structure, relative to the mass of an egg

"""
    init_substates_agent(p::AbstractParamCollection)
    
Initialize the agent substates, i.e. the state variables of a specific agent.
"""
function init_substates_agent(p::AbstractParamCollection)
    return ComponentArray( # initial states
        X_emb = Float64(p.agn.X_emb_int), # initial mass of vitellus
        S = Float64(p.agn.X_emb_int * X_EMB_INT_REL), # initial structure is a small fraction of initial reserve // mass of vitellus
        S_0 = Float64(p.agn.X_emb_int * X_EMB_INT_REL), # initial reference structure
        H = Float64(0), # maturity
        H_b = Float64(0), # maturity at birth (will be derived from model output)
        R = Float64(0), # reproduction buffer
        f_X = Float64(1), # scaled functional response 
        I_emb = Float64(0), # ingestion from vitellus
        I_p = Float64(0), # ingestion from external food resource
        I = Float64(0), # total ingestion
        A = Float64(0), # assimilation
        M = Float64(0), # somatic maintenance
        J = Float64(0), # maturity maintenance 
        Q = Float64(0), # cumulative dissipation flux

        D_G = MVector{length(p.spc.k_D_G), Float64}(zeros(length(p.spc.k_D_G))), # scaled damage | growth efficiency
        D_M = MVector{length(p.spc.k_D_M), Float64}(zeros(length(p.spc.k_D_M))), # scaled damage | maintenance costs 
        D_A = MVector{length(p.spc.k_D_A), Float64}(zeros(length(p.spc.k_D_A))), # scaled damage | assimilation efficiency
        D_R = MVector{length(p.spc.k_D_R), Float64}(zeros(length(p.spc.k_D_R))), # scaled damage | reproduction efficiency
        D_h = MVector{length(p.spc.k_D_h), Float64}(zeros(length(p.spc.k_D_h))), # scaled damage | hazard rate

        y_G = Float64(1.), # relative response | growth efficiency
        y_M = Float64(1.), # relative response | maintenance costs 
        y_A = Float64(1.), # relative response | assimilation efficiency
        y_R = Float64(1.), # relative response | reproduction efficiency
        h_z = Float64(0.), # hazard rate | chemical stressors
        h_S = Float64(0.),  # hazard rate | starvation

        cum_repro = Float64(1),
        age = Float64(0.),
        cohort = Float64(0),
        dead = Float64(0),
        causeofdeath = Float64(0)
    )
end

"""
    init_substates_global(p::AbstractParamCollection)
Initialize the global substates, i.e. the global state variables such as simulation time.
"""
function init_substates_global(p::AbstractParamCollection)
    return ComponentArray(
            X_p = Float64(p.glb.Xdot_in), # initial resource abundance equal to influx rate
            C_W = p.glb.C_W, # external stressor concentrations
            N_tot = p.glb.N0
        )
end

"""
    initialize_statevars(p::AbstractParamCollection, pagnt::ComponentVector{Float64})::ComponentArray 
For initialization of ODE simulator, initialize the component vector of state variables, `u`, based on common parameters `p`.
"""
function initialize_statevars(p::AbstractParamCollection)::ComponentArray 
    return ComponentArray(
        glb = init_substates_global(p),
        agn = init_substates_agent(p)
    )
end

