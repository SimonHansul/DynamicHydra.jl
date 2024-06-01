using Revise
using Hydra
default(palette = palette([:lightblue, :purple], 5), lw = 2)

p = DEBParamCollection()
p.glb.C_W = [0.]
p.spc

p.spc.k_D_R = [0.]
p.spc.k_D_M = [10.]
p.spc.e_M = [1.]
p.spc.b_M = [2.]

yhat = simulator(p)

using Plots, StatsPlots
@df yhat plot(:t, :S)

p.glb.C_W = [1.]
yhat = simulator(p)
@df yhat plot!(:t, :S)

