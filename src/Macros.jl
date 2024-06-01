function idcol!(Yhat::AbstractDataFrame, col::Symbol, val)
    Yhat[!,col] .= val
end

function idcol!(Yhat::Any, col::Symbol, val)
    for df in Yhat
        idcol!(df, col, val)
    end
end

concat_sims(Yhat::Vector{DataFrame}) = vcat(Yhat...)
concat_sims(Yhat::Vector{Any}) = ([vcat([Yhat[i][j] for i in eachindex(Yhat)]...) for j in eachindex(Yhat[1])]...,)

"""
    @replicates(simcall::Expr, nreps::Int64) 

Perform replicated runs of `simcall`, where `simcall` is a call to a simulator function. 

Example:

    spc = SpeciesParams(Z = Truncated(Normal(1, 0.1), 0, Inf)) # initialize default parameters with variable zoom factor
    yhat = @replicates Hydra.simulator(DEBParamCollection(spc = spc))) 10 # execute replicated runs to simulator

In this case, `yhat` will contain the output of 10 replicated simulations. For each replicate, the zoom factor is sampled from a truncated Normal distribution. 
`yhat` contains an additional column `replicate`. 
If yhat is a Vector{Any}, we assume that each element is an iterable of DataFrames.
"""
macro replicates(simcall::Expr, nreps::Int64)
    quote
        yhat = [] #DataFrame()

        for replicate in 1:$nreps
            yhat_i = $(esc(simcall))
            idcol!(yhat_i, :replicate, replicate)#yhat_i[!,:replicate] .= replicate
            push!(yhat, yhat_i)
        end
        yhat = concat_sims(yhat)
        yhat
    end

end
