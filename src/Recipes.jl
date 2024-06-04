
@userplot RugPlot
@recipe function f(h::RugPlot)
    if length(h.args)>1
        error("Rugplot does not support multiple arguments.")
    end
    @series begin
        seriestype := :scatter
        markershape := :vline
        x, y = h.args[1], repeat([0], length(h.args[1]))
    end
end

@userplot LinePlot
@recipe function f(h::LinePlot; estimator = mean)
    let x, y, q, mask
        
        if length(h.args)>=3
            x, y,  q = h.args
        else # if no quantile argument given at all => assume q=(0.05, 0.95)
            x, y, = h.args
            q = (.05, .95)
        end

        mask = isfinite.(x) .& isfinite.(y) 
        x = x[mask]
        y = y[mask]
        
        x_ag = []
        y_mn = []
        y_lo = []
        y_hi = []

        for (j,xi) in enumerate(unique(x))
            mask = x.==xi
            push!(x_ag, xi)
            push!(y_mn, estimator(y[mask]))
            push!(y_lo, quantile(y[mask], q[1]))
            push!(y_hi, quantile(y[mask], q[2]))
        end

        @series begin
            seriestype --> :path
            ribbon := @. (y_mn - y_lo, y_hi - y_mn)
            x_ag, y_mn
        end
    end
end

@userplot GroupedLinePlot
@recipe function f(h::GroupedLinePlot; estimator = mean)
    let x, y, group, q, mask
        # if no quantile argument given at all => assume q=(0.05, 0.95)
        if length(h.args)>=4
            x, y, group, q = h.args
        else
            x, y, group = h.args
            q = (.05, .95)
        end
        if typeof(group[1])<:Number
            mask = isfinite.(x) .& isfinite.(y) .& isfinite.(group) .& (ismissing.(x).==false) .& (ismissing.(y).==false) .& (ismissing.(group).==false)
        else
            mask = isfinite.(x) .& isfinite.(x) .& (ismissing.(x).==false) .& (ismissing.(y).==false)
        end
        x = x[mask]
        y = y[mask]
        group = group[mask]
        # aggregate values
        x_ag = []
        y_mn = []
        y_lo = []
        y_hi = []
        g_ag = []
        for (i,g) in enumerate(unique(group))
            for (j,xi) in enumerate(unique(x[group.==g]))
                mask = (group.==g) .& (x.==xi)
                push!(x_ag, xi)
                push!(y_mn, estimator(y[mask]))
                push!(y_lo, quantile(y[mask], q[1]))
                push!(y_hi, quantile(y[mask], q[2]))
                push!(g_ag, g)
            end
        end
        for (i,g) in enumerate(unique(g_ag))
            @series begin
                mask = g_ag .== g
                seriestype --> :path
                ribbon := (y_mn[mask] .- y_lo[mask], y_hi[mask] .- y_mn[mask])
                x_ag[mask], y_mn[mask]
            end
        end
    end
end

#=
Some additional useful functions...
=#

"""
Create shared y-labels for plot with grid layout, placing labels only on the left-most subplots.

"""
function gridylabel(label::String, numrows::Int64, numcols::Int64)::Matrix{String}
    labelstring = [label] # the label as vector
    numemptystrings = numcols-1 # the number of empty strings needed per row
    emptystrings = repeat([""], numemptystrings) # vector of empty strings
    labels = repeat(vcat(labelstring, emptystrings), numrows) # repeat for each row
    labels = hcat(labels...) # convert to matrix
    return labels
end

function gridxlabel(label::String, numrows::Int64, numcols::Int64)::Matrix{String}
    labelstring = [label] # the label as vector
    numemptyrows = numrows-1 # the number of rows with empty strings
    emptystrings = repeat([""], numemptyrows * numcols) # vector of empty strings
    lastrow = repeat(labelstring, numcols) # the last row, containing the actual labels
    labels = vcat(emptystrings, lastrow) # repeat for each row
    labels = hcat(labels...) # convert to matrix
    return labels
end

"""
Save plot plt with name, assuming known global TAG and plotsdir.

"""
function saveplt(plt, name)
    savefig(plt, plotsdir("$(TAG)_$(name).pdf"))
    savefig(plt, plotsdir("$(TAG)_$(name).png"))
end