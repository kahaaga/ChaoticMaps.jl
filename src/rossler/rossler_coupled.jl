using DifferentialEquations, Parameters

@with_kw struct RosslerCoupledParams
  a = 0.1 # parameter in the rossler attractor
  b = 0.1 # parameter in the rossler attractor
  c = 14.0 # parameter in the rossler attractor
  ϵ₁ = 0.0 # influence of subsystem 1 on subsystem 2
  ϵ₂ = 0.0 # influence of subsystem 2 on subsystem 1
  ω₁ = 1+0.015 # the frequency of the first system
  ω₂ = 1-0.015 # the frequency of the second system
end

"""
    rossler_coupled_symbolic! = @ode_def RosslerCoupled begin
        # First Rossler system
        dx₁ = ω₁*(-y₁) - z₁ + ϵ₂*(x₁ - x₂)
        dy₁ = ω₁*x₁ + a*y₁
        dz₁ = b + z₁*(x₁ - c)

        dx₂ = ω₂*(-y₂) - z₂ + ϵ₁*(x₂ - x₁)
        dy₂ = ω₂*x₂ + a*y₂
        dz₂ = b + z₂*(x₂ - c)
    end a b c ϵ₁ ϵ₂ ω₁ ω₂

Symbolic representation of a coupled Rossler system. The system consists of two
separate subsystems, each being a Rossler attractor. The subsystems influence
each other through variables x₁ and x₂, with ϵ₁` and `ϵ₂` controlling the
relative influences.

# Examples

```julia-repl
# Define a set of Rossler parameters. We're using the default here, but you can
# also assign specific parameters as RosslerCoupledParams(ϵ₂ = 0.1, ω₁ = 1.1);
# for the remaining parameters, default values are used.
rp = RosslerCoupledParams()
# Represent params as tuple (needed for the symbolically defined system)
params = [getfield(rp, x) for x in fieldnames(rp)]

# Set up problem and solve it
u0 = rand(6) # initial conditions
tspan = (0.0, 300.0)
prob = ODEProblem(rossler_coupled_symbolic!, rand(6), tspan, p)
sol = solve(prob, Rosenbrock23(), saveat = 0.05)

# Plot stuff
using Plots; plotlyjs()

# 3D plot
# system 1
p1 = plot(sol, vars=(1,2,3), width = 0.8, color = "black", denseplot = false)
# system 2
plot!(p1, sol, vars=(4,5,6), width = 0.8, color = "red", denseplot = false)

# Plot x₁ and x₂ (note: var 0 is time)
p2 = plot(sol, vars=(0, 1), width = 0.8, color = "black", denseplot = false)  # x₁
plot!(p2, sol, vars=(0, 4), width = 0.8, color = "red", denseplot = false)   # x₂
```
"""
rossler_coupled_symbolic! = @ode_def RosslerCoupledSystem begin
    # First Rossler system
    dx₁ = ω₁*(-y₁) - z₁ + ϵ₂*(x₁ - x₂)
    dy₁ = ω₁*x₁ + a*y₁
    dz₁ = b + z₁*(x₁ - c)

    dx₂ = ω₂*(-y₂) - z₂ + ϵ₁*(x₂ - x₁)
    dy₂ = ω₂*x₂ + a*y₂
    dz₂ = b + z₂*(x₂ - c)
end a b c ϵ₁ ϵ₂ ω₁ ω₂

"""
    solve_rossler_coupled(; rp = RosslerCoupledParams(),
                            tend = 300.0,
                            timestep = 0.5,
                            u0 = rand(6))

Solves the coupled rossler system (for more details, see docs for
`rossler_coupled_symbolic!``).

The solution goes from time zero until `tend` and is sampled at regular time
intervals of width `timestep`.

`u0` is a 6-element vector of initial conditions (3 for each subsystem), and
`rp` is an instance of `RosslerCoupledParams`.
"""
function solve_rossler_coupled(; rp = RosslerCoupledParams(),
                                tend = 300.0,
                                timestep = 0.5,
                                u0 = rand(6))

  # Define a set of Rossler parameters. We're using the default here, but you can
  # also assign specific parameters as RosslerCoupledParams(ϵ₂ = 0.1, ω₁ = 1.1);
  # for the remaining parameters, default values are used.
  rp = RosslerCoupledParams()

  # Represent params as tuple (needed for the symbolically defined system)
  params = [getfield(rp, x) for x in fieldnames(rp)]

  # Set up problem and solve it
  tspan = (0.0, tend)
  prob = ODEProblem(rossler_coupled_symbolic!, u0, tspan, params)

  solve(prob, Rosenbrock23(), saveat = timestep)
end
