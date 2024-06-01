
@enum ReturnType dataframe odesolution matrix # possible return types

function abstractsimulator(
    p::AbstractParamCollection,
    model, 
    AgentParamType::DataType;
    alg = Tsit5(),
    returntype::ReturnType = dataframe,
    kwargs...
    )::Union{DataFrame,ODESolution}

    p.agn = AgentParamType(p.spc) # initialize agent parameters incl. individual variability
    
    u = initialize_statevars(p)
    prob = ODEProblem(model, u, (0, p.glb.t_max), p) # define the problem
    sol = solve(prob, alg; kwargs...) # get solution to the IVP

    if returntype == dataframe
        return sol_to_df(sol) # convert solution to dataframe
    end
    
    if returntype == matrix
        return sol_to_mat(sol) # convert solution to matrix
    end

    if returntype == odesolution
        return sol # directly return the ODESolution object
    end

end

"""
    simulator(
        p::Ref{DEBParamCollection}; 
        alg = Tsit5(),
        saveat = 1, 
        reltol = 1e-6,
        kwargs...
        )::DataFrame

Run the DEBBase model from a `DEBParamCollection` instance. <br>
"""
function simulator(
    p::DEBParamCollection; 
    model = DEBODE!,
    AgentParamType::DataType = AgentParams,
    kwargs...
    )

    abstractsimulator(
        p, 
        model,
        AgentParamType;     
        alg = Tsit5(),
        saveat = 1, 
        reltol = 1e-6,
        kwargs...
        )
end

