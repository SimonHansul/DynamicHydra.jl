# Solvers.jl
# Currently this only contains the Euler scheme. 
# We could think about more sophisticated integration between the ABM and DifferentialEquations.jl, as described in https://juliadynamics.github.io/jl/v5.7/examples/diffeq/. 
# Currently not clear however how to best combine this with stochastic processes in the ABM. TBD.

"""
    euler(du::Float64, u::Float64, dt::Float64)::Float64

Apply Euler scheme.
"""
function euler(du::Float64, u::Float64, dt::Float64)::Float64
    return u + du * dt
end


"""
    defeuler(dt::Float64)

Define an Euler function for fixed `dt`. Used by `ABM`.
"""
function defeuler(dt::Float64)
    function euler_dt(du::Float64, u::Float64)
        return euler(du, u, dt)
    end
    return euler_dt
end