function add_idcol(Yhat::AbstractDataFrame, col::Symbol, val)
    Yhat[!,col] .= val
    return Yhat
end

concat_sims(Yhat::Vector{DataFrame}) = vcat(Yhat...)

"""
    @replicates(simcall::Expr, nreps::Int64) 

Perform replicated runs of `simcall`, where `simcall` is a call to a simulator function. 

Example:

    spc = SpeciesParams(Z = Truncated(Normal(1, 0.1), 0, Inf)) # initialize default parameters with variable zoom factor
    yhat = @replicates Hydra.simulator(ODEParamCollection(spc = spc))) 10 # execute replicated runs to simulator

In this case, `yhat` will contain the output of 10 replicated simulations. For each replicate, the zoom factor is sampled from a truncated Normal distribution. 
`yhat` contains an additional column `replicate`. 
If yhat is a Vector{Any}, we assume that each element is an iterable of DataFrames.
"""
macro replicates(simcall::Expr, nreps::Int64)
    quote
        yhat = DataFrame[]

        for replicate in 1:$nreps
            yhat_i = $(esc(simcall))
            yhat_i = add_idcol(yhat_i, :replicate, replicate)#yhat_i[!,:replicate] .= replicate
            push!(yhat, yhat_i)
        end
        yhat = concat_sims(yhat)
        yhat
    end

end
