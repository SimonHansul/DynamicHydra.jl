if occursin("terminal", abspath(PROGRAM_FILE))
    @info("Loading packages")
    using Pkg; Pkg.activate("test")
    using Plots, StatsPlots, Plots.Measures
    using SHUtils
    using DataFrames, DataFramesMeta
    using ProgressMeter
    default(leg = false, lw = 1.5)
    using DEBFigures
    using OrdinaryDiffEq

    using Revise
    using DynamicHydra
end
TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")

pmoa = "M"
p = DEBParamCollection()
p.glb.t_max = 42.
p.spc.k_D_G = [0.]
p.spc.k_D_M = [10.]
p.spc.k_D_A = [0.]
p.spc.k_D_R = [0.]
p.spc.k_D_h = [0.]

p.spc.e_M = [1.]
p.spc.b_M = [2.]

#=
Simulate single stressors with different PMoAs
=#
# FIXME: negative damage...
# does not seem related to S at all
#
using DynamicHydra
begin   
    out = DataFrame()
    pmoas = ["G", "M", "A", "R"]
    for pmoa in pmoas
        for C_W in round.(10 .^ range(log10(0.1), log10(1.), length = 5), sigdigits = 2)
            glb = GlobalParams(t_max = 42., C_W = [C_W])
            spc = SpeciesParams(
                kappa = 0.538,
                k_D_G = [0.5], 
                k_D_M = [0.5], 
                k_D_A = [0.5], 
                k_D_R = [0.5], 
                k_D_h = [0.5], 
                e_G = [1.],
                e_M = [1.],
                e_A = [1.],
                e_R = [1.],
                e_h = [1.],
                b_G = [2.],
                b_M = [2.],
                b_A = [2.],
                b_R = [2.],
                b_h = [2.])
            p = DEBParamCollection(glb = glb, spc = spc)
            isolate_pmoas!(p.spc, [pmoa])
            out_zj = simulator(p)
            out_zj[!,:C_W] .= C_W
            out_zj[!,:pmoa] .= pmoa
            append!(out, out_zj)
        end
    end
    out = DynamicHydra.relative_response(out, [:S, :R], :C_W; groupby_vars = [:t, :pmoa])

    plt = plot(
        layout = (2,4), title = hcat(pmoas...), 
        ylim = (0, 1.01),
        xlabel = "t", ylabel = hcat([gridylabel("y_S", 1, 4), gridylabel("y_R", 1, 4)]...), 
        size = (800,450), 
        bottommargin = 5mm, leftmargin = 2.5mm, lw = 2,
        leg = [false false false true false false false false],
        legtitle = "C_W / e", legendtitlefontsize = 10
        )

    for (j,pmoa) in enumerate(pmoas)
        @df @subset(out, :pmoa .== pmoa) groupedlineplot!(plt, :t, :y_S, :C_W, subplot = j, label = hcat(unique(:C_W)...))
        @df @subset(out, :pmoa .== pmoa) groupedlineplot!(plt, :t, :y_R, :C_W, subplot = 4+j, label = hcat(unique(:C_W)...))
    end

    display(plt)
    savefig(plt, "plots/$(TAG).png")
end

