
#=
# Simulator function for ODE-based DEB model
=#

#@enum ReturnType dataframe odesolution matrix # possible return types

#function abstractsimulator(
#    p::AbstractParamCollection,
#    model, 
#    AgentParamType::DataType;
#    alg = Tsit5(),
#    returntype::ReturnType = dataframe,
#    kwargs...
#    )::Union{DataFrame,ODESolution}
#
#    p.agn = AgentParamType(p.spc) # initialize agent parameters incl. individual variability
#    
#    u = initialize_statevars(p)
#    prob = ODEProblem(model, u, (0, p.glb.t_max), p) # define the problem
#    sol = solve(prob, alg; kwargs...) # get solution to the IVP
#
#    if returntype == dataframe
#        return sol_to_df(sol) # convert solution to dataframe
#    end
#    
#    if returntype == matrix
#        return sol_to_mat(sol) # convert solution to matrix
#    end
#
#    if returntype == odesolution
#        return sol # directly return the ODESolution object
#    end
#
#end

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

