"""
    calc_SL_max(spc::AbstractSpeciesParams)::Float64
Calculate maximum structural length slmax [m^(1/3)]
"""
function calc_SL_max(spc::AbstractSpeciesParams)::Float64
    return ((spc.kappa * spc.Idot_max_rel * spc.eta_IA) / spc.k_M)
end

"""
    calc_S_max(spc::AbstractParams)::Float64
Calculate maximum structural mass smax [m]
"""
function calc_S_max(spc::AbstractSpeciesParams)::Float64
    return calc_SL_max(spc)^3
end

"""
    k_J!(spc::AbstractSpeciesParams)::Nothing
Set the maturity maintenance rate constant, 
assuming that the cumulative investment into maturity maintenance 
equals the cumulative investment into somatic maintenance.
"""
function k_J!(spc::AbstractSpeciesParams)::Nothing
    spc.k_J = ((1 - spc.kappa) / spc.kappa) * spc.k_M

    return nothing
end