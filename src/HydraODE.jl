"""
    HydraODE!(du, u, p, t)
Definition of  a Hydra model as a system of ordinary differential equations. 
"""
function HydraODE!(du, u, p, t)::Nothing

    for func! in p.glb.odefuncs # calculate the global derivatives
        func!(du, u, p, t) 
    end

    for func! in p.spc.odefuncs # calculate the species-specific derivatives
        func!(du, u, p, t)
    end

    return nothing
end

"""
    simulator(
        theta::AbstractParamCollection; 
        AgentParamType::DataType = AgentParams,
        kwargs...
        )::DataFrame

Run a from any `DEBParamCollection` instance. <br>
The `AgentParamType` only has to be changed if a new kind of agent variability is introducted. <br>
kwargs are passed down to `OrdinaryDiffEq.solve()`. <br>
"""
function simulator(
    p::AbstractParamCollection; 
    alg = Tsit5(),
    saveat = 1.,
    AgentParamType::DataType = AgentParams,
    kwargs...
    )

    p.agn = AgentParamType(p.spc) # initialize agent parameters incl. individual variability
    u = initialize_statevars(p)
    prob = ODEProblem(HydraODE!, u, (0, p.glb.t_max), p) # define the initial value problem
    sol = solve(prob, alg; saveat = saveat, kwargs...) # get solution to the IVP
    return sol_to_df(sol)
end

