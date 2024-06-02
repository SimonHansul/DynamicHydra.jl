"""
Clip negative values.
"""
function clipneg(x::Float64)
    return sig(x, 0., 0., x)
end

"""
Sigmoid switch function. 
`y_left` and `y_right` are the function values left and right of the threshold `x_thr`.

"""
@inline function sig(
    x::Float64, 
    x_thr::Float64,
    y_left::Float64, 
    y_right::Float64; 
    beta::Float64 = 1e16
    )::Float64

    return 1 / (1 + exp(-beta*(x - x_thr))) * (y_right - y_left) + y_left
end

@inline function functional_response(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Float64

    let X_V = u_glb.X_p / p_glb.V_patch # convert food abundance to concentration
        return X_V / (X_V + p_agn.K_X) # calculate type II functional response
    end
end

"""
Calculate ingestion rate. 
Embryos (X_emb <= 0) take up resources from the vitellus X_emb. 
Juveniles and adults (X_emb > 0) feed on the external resource X_pcmn.

"""
@inline function Idot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing

    du_agn.f_X = functional_response(du_glb, du_agn, u_glb, u_agn, p_glb, p_spc, p_agn,t)
    
    du_agn.I_emb = sig( # uptake from vitellus
        du_agn.X_emb, # uptake from vitellus depends on mass of vitellus
        0., # the switch occurs when vitellus is used up 
        0., # when the vitellus is used up, there is no uptake
        (Complex(du_agn.S)^(2/3)).re * p_agn.Idot_max_rel; # when the vitellus is not used up, uptake from vitellus occurs
        beta = 1e20 # for switches around 0, we need very high beta values
        )

    du_agn.I_p = sig( # uptake from external resource p
        du_agn.X_emb, # ingestion from external resource depends on mass of vitellus
        0., # the switch occurs when the vitellus is used up  
        du_agn.f_X * p_agn.Idot_max_rel * (Complex(du_agn.S)^(2/3)).re, # when the vitellus is used up, ingestion from the external resource occurs
        0.; # while there is still vitellus left, there is no uptake from the external resource
        beta = 1e20 # again we have a switch around 0, requiring very high beta
        )

    du_agn.I = du_agn.I_emb + du_agn.I_p

    du_glb.X_p -= du_agn.I_p
    du_agn.X_emb = -du_agn.I_emb

    return nothing
end

"""
Assimilation flux
"""
@inline function Adot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real)::Nothing

    du_agn.A = du_agn.I * p_spc.eta_IA * du_agn.y_A

    return nothing
end

"""
Somatic maintenance flux
"""
@inline function Mdot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real)::Nothing

    du_agn.M = du_agn.S * p_spc.k_M * du_agn.y_M

    return nothing
end

"""
Maturity maintenance flux

"""
@inline function Jdot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing

    du_agn.J = du_agn.H * p_spc.k_J * du_agn.y_M

    return nothing
end

"""
Positive somatic growth
"""
@inline function Sdot_positive(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Float64

    return p_spc.eta_AS * du_agn.y_G * (p_spc.kappa * du_agn.A - du_agn.M)
end

"""
Negative somatic growth
"""
@inline function Sdot_negative(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Float64 

    return -(du_agn.M / p_spc.eta_SA - p_spc.kappa * du_agn.A)
end

function Sdot(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Float64 

    return sig(
        p_spc.kappa * du_agn.A, # growth depends on maintenance coverage
        du_agn.M, # switch occurs based on maintenance costs
        Sdot_negative(du_glb, du_agn, u_glb, u_agn, p_glb, p_spc, p_agn,t), # left of the threshold == maintenance costs cannot be covered == negative growth
        Sdot_positive(du_glb, du_agn, u_glb, u_agn, p_glb, p_spc, p_agn, t)# right of the threshold == maintenance costs can be covered == positive growth
    )
end

"""
Somatic growth, including application of shrinking equation.
"""
@inline function Sdot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing

    du_agn.S = Sdot(du_glb, du_agn, u_glb, u_agn, p_glb, p_spc, p_agn,t)

    return nothing
end


"""
Maturation flux. 

Maturity is dissipated energy and can therefore not be burned to cover maintenance costs. 

Currently there are no consequences for an organism not covering maturity maintenance. 
If this turns out to be an issue, we might consider to add a damage pool ``D_H``,
where the amount of maturity maintenance that could not be covered is accumulated. 
This might then lead to a fitness penalty depending on ``D_H``, for example in the form of additional 
mortality or embryonic hazard. 

Such rules are however likely species-specific and should be evaluated in the light of a more precise problem definition.
"""
@inline function Hdot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing 

    du_agn.H = sig(
        du_agn.H, # maturation depends on maturity
        p_agn.H_p, # switch occurs at maturity at puberty H_p
        clipneg(((1 - p_spc.kappa) * du_agn.A) - du_agn.J), # maturation for embryos and juveniles
        0., # maturation for adults
    )

    return nothing
end

"""
Update the current estimate of H_b. 
The current estimate of maturity at birth is equal to current maturity for embryos, 
and will be fixed to the current value upon completion of embryonic development.
This way, we obtain H_b as an implied trait and can use it later (for example in `abj()`).
"""
@inline function H_bdot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing

    du_agn.H_b = sig(
        du_agn.X_emb, # estimate depends on embryonic state
        0., # switch occurs when vitellus is gone
        0., # post-embryonic: H_b stays fixed
        du_agn.H # embryonic: H_b tracks H
    )

    return nothing
end

"""
Reproduction flux.
"""
@inline function Rdot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing 

    du_agn.R = sig(
        du_agn.H, # reproduction depends on maturity
        p_agn.H_p, # switch occurs at maturity at puberty H_p
        0., # reproduction for embryos and juveniles
        clipneg(du_agn.y_R * p_spc.eta_AR * ((1 - p_spc.kappa) * du_agn.A - du_agn.J)) # reproduction for adults
    )

    return nothing
end

"""
    Qdot!(
        du::ComponentArray,
        u::ComponentArray,
        p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
        t::Real
        )::Nothing 

Calculation of the total dissipation flux, equal to the sum of maintenance costs and overheads paid for assimilation, mobilization, maturation, growth and reproduction.
"""
function Qdot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing 

    # dissipation fluxes for the individual processes
    let Qdot_A, Qdot_S, Qdot_C, Qdot_R
        
        Qdot_A = du.I * (1 - p_spc.eta_IA)
        Qdot_S = du.S >= 0 ? du.S * (1 - p_spc.eta_AS) / p_spc.eta_AS : du.S * (p_spc.eta_SA - 1)
        Qdot_R = du.R * (1 - p_spc.eta_AR) / p_spc.eta_AR
        
        du.Q =  Qdot_A + Qdot_S + Qdot_C + Qdot_R + du.M + du.J + du.H
    end

    return nothing
end

"""
    C_Wdot!(
        du::ComponentVector,
        u::ComponentVector,
        p_glb::AbstractGlobalParams,
p_spc::AbstractSpeciesParams,
p_agn::AbstractAgentParams,
        t::Real
        )::Nothing 

Change in external concentrations. 
Currently simply returns zeros because time-variable exposure is not yet implemented.
"""
@inline function C_Wdot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing 

    du_glb.C_W = zeros(length(u_glb.C_W)) # constant exposure : derivative is 0

    return nothing
end

"""
Change in environmental resource abundance for simulation of a chemostat.
"""
function X_pdot_chemstat!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing 

    du.X_p = p_glb.Xdot_in - p_glb.k_V * u.X_p

    return nothing
end

"""
Constant concentration of external chemical stressor:
"""
function C_Wdot_const!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing

    du.C_W = zeros(length(u.C_W))

    return nothing
end


"""
TK for spc-TKTD model, including effect of surface area to volume ratio and dilution by growth. 
If `D` and is given as a Vector, TK is only stressor-specific but not PMoA-specific. 
"""
@inline function Ddot!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing 

    for z in eachindex(u_glb.C_W)
        @inbounds du_agn.D_G[z] = sig(du_agn.X_emb, 0., p_spc.k_D_G[z] * (u_glb.C_W[z] - du_agn.D_G[z]), 0.)
        @inbounds du_agn.D_M[z] = sig(du_agn.X_emb, 0., p_spc.k_D_M[z] * (u_glb.C_W[z] - du_agn.D_M[z]), 0.)
        @inbounds du_agn.D_A[z] = sig(du_agn.X_emb, 0., p_spc.k_D_A[z] * (u_glb.C_W[z] - du_agn.D_A[z]), 0.)
        @inbounds du_agn.D_R[z] = sig(du_agn.X_emb, 0., p_spc.k_D_R[z] * (u_glb.C_W[z] - du_agn.D_R[z]), 0.)
        @inbounds du_agn.D_h[z] = sig(du_agn.X_emb, 0., p_spc.k_D_h[z] * (u_glb.C_W[z] - du_agn.D_h[z]), 0.)
        
        #@inbounds du.D_G[z] = sig(u.X_emb, 0., p_spc.k_D_G[z] * (calc_SL_max(p.spc) / (Complex(u.S)^(1/3)).re) * (u.C_W[z] - u.D_G[z]) - u.D_G[z] * (du.S / u.S), 0.)
        #@inbounds du.D_M[z] = sig(u.X_emb, 0., p_spc.k_D_M[z] * (calc_SL_max(p.spc) / (Complex(u.S)^(1/3)).re) * (u.C_W[z] - u.D_M[z]) - u.D_M[z] * (du.S / u.S), 0.)
        #@inbounds du.D_A[z] = sig(u.X_emb, 0., p_spc.k_D_A[z] * (calc_SL_max(p.spc) / (Complex(u.S)^(1/3)).re) * (u.C_W[z] - u.D_A[z]) - u.D_A[z] * (du.S / u.S), 0.)
        #@inbounds du.D_R[z] = sig(u.X_emb, 0., p_spc.k_D_R[z] * (calc_SL_max(p.spc) / (Complex(u.S)^(1/3)).re) * (u.C_W[z] - u.D_R[z]) - u.D_R[z] * (du.S / u.S), 0.)
        #
        #@inbounds du.D_h[z] = sig(u.X_emb, 0., p_spc.k_D_h[z] * (calc_SL_max(p.spc) / (Complex(u.S)^(1/3)).re) * (u.C_W[z] - u.D_h[z]) - u.D_h[z] * (du.S / u.S), 0.)
    end

    return nothing
end

@inline function age!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing 

    du_agn.age = 1.

    return nothing
end

"""
Response to chemical stressors, assuming independent action for mixtures.
"""
@inline function y_z!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing 

    @inbounds du_agn.y_G = prod([p_spc.drc_functs_G[z](du_agn.D_G[z], (p_spc.e_G[z], p_spc.b_G[z])) for z in 1:length(u_glb.C_W)]) # combined relative responses for sublethal effects per PMoA
    @inbounds du_agn.y_M = prod([p_spc.drc_functs_M[z](du_agn.D_M[z], (p_spc.e_M[z], p_spc.b_M[z])) for z in 1:length(u_glb.C_W)])
    @inbounds du_agn.y_A = prod([p_spc.drc_functs_A[z](du_agn.D_A[z], (p_spc.e_A[z], p_spc.b_A[z])) for z in 1:length(u_glb.C_W)])
    @inbounds du_agn.y_R = prod([p_spc.drc_functs_R[z](du_agn.D_R[z], (p_spc.e_R[z], p_spc.b_R[z])) for z in 1:length(u_glb.C_W)])

    @inbounds du_agn.h_z = sum([p_spc.drc_functs_h[z](du_agn.D_h[z], (p_spc.e_h[z], p_spc.b_h[z])) for z in 1:length(u_glb.C_W)]) # hazard rate

    return nothing
end


"""
Hazard rate under starvation
"""
@inline function h_S!(
    du_glb::ComponentArray,
    du_agn::ComponentArray,
    u_glb::ComponentArray,
    u_agn::ComponentArray,
    p_glb::AbstractGlobalParams,
    p_spc::AbstractSpeciesParams,
    p_agn::AbstractAgentParams,
    t::Real
    )::Nothing 

    du_agn.h_S = LL2h(du_agn.S / du_agn.S_0, (p_spc.e_S, -p_spc.b_S))

    return nothing
end

"""
Maturity-driven metabolic acceleration from birth to maturity threshold `H_j` (metamorphosis). 
We assume that some baseline parameter `p` has value `p_b` at birth and `p_j` at metamorphosis.
Between birth and metamorphosis, the current value of `p` is the maturity-weighted mean of `p_b` and `p_j`.
"""
function Hbj(H::Float64, X_emb::Float64, H_b::Float64, H_j::Float64, p_b::Float64, p_j::Float64)::Float64
    w_b = (H_j - H) / (H_j - H_b) # weight for p_b
    w_j = 1 - w_b # weight for p_j
    p_bj = mean([p_b, p_j], Weights([w_b, w_j])) # p_bj, i.e. value between birth and maturity
    
    p = Hydra.sig( # post-metamorphosis: value stays constant at p_j
        H,
        H_j,
        Hydra.sig( # embryonic: value stays constant at p_b
            X_emb, 
            0., 
            p_bj, 
            p_b),
        p_j
    )

    return p
end


"""
    reproduce_opportunistic!(agent::AbstractAgent, abm::AbstractABM)::Nothing

Reproduction according to an opportunistic reproduction strategy. 
Offspring is created whenever there is enough energy in the reproduction buffer.
"""
function reproduce_opportunistic!(agent::AbstractAgent, abm::AbstractABM)::Nothing
    
    let num_offspring = trunc(agent.du_agn.R / agent.p_agn.X_emb_int) # calculate the number of offspring produced
        for _ in 1:num_offspring # for the given number of offspring agents
            offspring = abm.p_glb.AgentType(abm) # initialize a new agent
            offspring.du_agn.cohort = agent.du_agn.cohort + 1 # the offspring agent is part of a new cohort
            push!(abm.agents, offspring) # add offspring to agents vector
            agent.du_agn.R -= agent.p_agn.X_emb_int # remove energy from reproduction buffer
        end
        agent.du_agn.cum_repro += num_offspring
    end

    return nothing
end

"""
    ingestion!(agent::AbstractAgent, abm::AbstractABMagent)::Nothing

Update resource abundance in model `abm` based on ingestion flux of `agent`.
"""
function ingestion!(agent::AbstractAgent, abm::AbstractABM)::Nothing
    
    abm.u.X_p = max(0, abm.u.X_p - agent.du_agn.I_p * abm.dt)
    
    return nothing
end

"""
    stochastic_death!(agent::AbstractAgent, abm::AbstractABM, h_x::Float64)::Nothing

Apply stochastic death model to agent, with respect to hazard rate `h_x`. 

Records the cause of death encoded in a number, given by keyword argument `causeofdeath`.
"""
function stochastic_death!(agent::AbstractAgent, abm::AbstractABM, h_x::Float64; causeofdeath = 0.)::Nothing
    if rand() >= exp(-h_x / abm.dt)
        agent.du_agn.dead = 1.
        agent.du_agn.causeofdeath = causeofdeath
    end

    return nothing
end

function die!(agent::AbstractAgent, abm::AbstractABM)::Nothing
    stochastic_death!(agent, abm, agent.du_agn.h_S; causeofdeath = 1.) # starvation mortality from loss of structure
    return nothing
end

function N_tot!(abm::AbstractABM)::Nothing
    abm.u.N_tot = length(abm.agents)
    return nothing
end

"""
    HydraODE!(du, u, u, p, t)
Generic definition of a Hydra ODE system. 
"""
function HydraODE!(du, u, p, t)::Nothing

    @unpack p_glb, p_spc, p_agn = p
    @unpack u_glb, u_agn = u
    @unpack du_glb, du_agn = du

    for odefunc! in p_glb.odefuncs
        odefunc!(du_glb, du_agn, u_glb, u_agn, p_glb, p_spc, p_agn, t)
    end

    for odefunc! in p_spc.odefuncs 
        odefunc!(du_glb, du_agn, u_glb, u_agn, p_glb, p_spc, p_agn, t)
    end
    
    return nothing
end