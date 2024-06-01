

"""
Definition of basic ABM object.
Currently assumes that all agents are represented by a single agent type `AgentType`.
"""
mutable struct ABM <: AbstractABM
    p::AbstractParamCollection # parameter collection
    t::Float64 # current simulation time
    dt::Float64 # timestep
    euler::Function # definition of the euler function for the given dt
    saveat::Float64
    
    u::ComponentVector # global state variables
    du::ComponentVector # global derivatives

    AgentID_count::Int64 # cumulative number of agents in the simulation
    agents::AbstractVector # population of agents
    aout::AbstractVector # agent output
    mout::AbstractVector # model output

    """
    Instantiate ABM from param collection `p`. 
    """
    function ABM(
        p::A; 
        dt = 1/24, 
        saveat = 1,
        execute = true
        ) where A <: AbstractParamCollection

        abm = new() # initialize ABM object

        abm.p = p # store the parameters
        abm.t = 0. # initialize simulation time
        abm.dt = dt # timestep is given as keyword argument
        abm.euler = defeuler(dt) # Euler function for the given dt
        abm.saveat = saveat
        
        abm.u = init_substates_global(p) # initialize the global substates
        abm.du = similar(abm.u) # initialize the global derivatives
        abm.du.X_p = 0.
        abm.du.C_W .= 0.
        
        abm.AgentID_count = 0 # set agent count
        initialize_agents!(abm) # initailize the population of agents
        abm.aout = [] # initialize agent output
        abm.mout = [] # initialize model output

        if execute
            run!(abm) # run the model

            return prepare_output(abm) # prepare output in dataframes and return
        else
            return abm
        end
    end
end



"""
The definition of a generic DEBAgent.
"""
mutable struct DEBAgent <: AbstractAgent
    p::AbstractParamCollection
    u::ComponentVector
    du::ComponentVector
    AgentID::Int64

    """
        DEBAgent(abm::ABM)
    
    Initialization of a generic DEB agent within ABM structure.
    """
    function DEBAgent(abm::ABM)
        a = new()

        a.p = copy(abm.p) # parameters; TODO: #29 avoid copying everything
        a.p.agn = AgentParams(a.p.spc) # assign agent-level parameters (induces individual variability)
        initialize_statevars!(a, abm) # initialize agent-level state variables
        a.du = similar(a.u) # initialize agent-level derivatives
        a.du.glb = abm.du # reference to global derivatives
        a.du.agn = copy(a.u.agn) # derivatives of agent substates have same shape as states
        a.AgentID = abm.AgentID_count # assign AgentID
        abm.AgentID_count += 1 # increment AgentID counter
        
        return a
    end
end

"""
    init_agents!(abm::AbstractABM)::nothing
Initialize the population of 
"""
function initialize_agents!(abm::AbstractABM)::Nothing

    abm.agents = [] # initialize a vector of agents with undefined values and defined length

    for i in 1:abm.p.glb.N0 # for the number of initial agents
        push!(abm.agents, DEBAgent(abm)) # initialize an agent and add it to the vector of agents
    end

    return nothing
end


"""
    prepare_output(abm::AbstractABM)

Prepare ABM output, converting Vectors of ComponentArrays into DataFrames. 
Two separate DataFrames are returned containing global and agent-level state variables, respectively.
"""
function prepare_output(abm::AbstractABM)

    mout = DataFrame(hcat([[x.t] for x in abm.mout]...)', [:t]) |> 
    x-> hcat(x, DataFrame(hcat([m.u for m in abm.mout]...)', extract_colnames(abm.mout[1].u)))
    
    if abm.p.glb.recordagentvars
        aout = DataFrame(hcat([[x.t, x.AgentID] for x in abm.aout]...)', [:t, :AgentID]) |> 
        x-> hcat(x, DataFrame(hcat([a.u for a in abm.aout]...)', extract_colnames(abm.aout[1].u)))
    else
        aout = DataFrame()
    end

    return mout, aout
end

function update_glb!(agent::AbstractAgent, abm::AbstractABM)
    map!(x -> x, agent.u.glb, abm.u) 
end

"""
step!(agent::AbstractAgent, abm::AbstractABM; odefuncs::Vector{Function}, rulefuncs::Vector{Function})

Definition of agent step. \n
This definition is generic, so that the function body does not have to be modified 
in order to simulate different models.
"""
function step!(agent::DEBAgent, abm::ABM)
    du, u, p = agent.du, agent.u, agent.p # unpack agent substates, derivatives and parameters
    t = abm.t

    update_glb!(agent, abm) # update references to global state variables
    
    for func! in agent.p.spc.rulefuncs # apply all rule-based functions
        func!(agent, abm)
    end

    for func! in agent.p.spc.odefuncs # apply all ODE-based functions
        func!(du, u, p, t) 
    end


    map!(abm.euler, u.agn, du.agn, u.agn) # apply the euler scheme to agent substates
end

"""
record!(agent::DEBAgent, abm::ABM)

Record agent output (only if `p.glb.recordagentvars == true`).
"""
function record!(agent::DEBAgent, abm::ABM)
    if abm.p.glb.recordagentvars && isapprox(abm.t%abm.saveat, 0, atol = abm.dt)
        push!(abm.aout, (t = abm.t, AgentID = agent.AgentID, u = copy(agent.u.agn)))
    end
end


"""
step!(
    abm::ABM; 
    )

Execution of a generic ABM step, following the schedule: 

1. Shuffle agents
2. Calculate global derivatives
3. Calculate agent derivatives
4. Update agent states
5. Record agent states
6. Update global states
7. Record global states

It is important that the agent steps occur between calculation of the global derivatives and 
updating the global states, because global derivatives which may be influenced by agent derivatives 
are initialized during calculation of the global states and mutated by the  
Changing this order will lead to erroneous computation of the interaction between agents and the environment.
"""
function step!(abm::ABM)

    du, u, p = abm.du, abm.u, abm.p
    t = abm.t

    shuffle!(abm.agents)

    for func! in abm.p.glb.odefuncs # execute global ODE-based functions
        func!(du, u, p, t)
    end

    for func! in abm.p.glb.rulefuncs # execute global rule-based functions
        func!(abm)
    end

    for a in abm.agents # for every agent
        step!(a, abm) # execute agent steps
        record!(a, abm) # record agent data
    end

    map!(abm.euler, u, du, u) # apply the euler scheme to global states
    record!(abm) # record global states
    filter!(a -> a.u.agn.dead == false, abm.agents) # remove dead agents

    abm.t += abm.dt
end

"""
record!(abm::ABM)

Record model-level output.
"""
function record!(abm::ABM)
    if isapprox(abm.t%abm.saveat, 0, atol = abm.dt)
        push!(abm.mout, (t = abm.t, u = copy(abm.u)))
    end
end




"""
Run the ABM.
"""
function run!(abm::AbstractABM)
    while abm.t <= abm.p.glb.t_max
        step!(abm)
    end
end
