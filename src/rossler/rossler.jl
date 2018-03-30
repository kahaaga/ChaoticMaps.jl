using DifferentialEquations, Parameters

"""
    rossler_symbolic! = @ode_def Rossler begin
        # First Rossler system
        dx = -y - z
        dy = x + a*y
        dz = b + z*(x - c)
    end a b c

Symbolic representation the Rossler system.


# Examples
```julia-repl

# Define a set of Rossler parameters. We're using the default here, but you can
# also assign specific parameters as RosslerCoupledParams(ϵ₂ = 0.1, ω₁ = 1.1);
# for the remaining parameters, default values are used. Remember to represent
# params as tuple, which is needed for the symbolically defined system.

rp = RosslerParams()
params = [getfield(rp, x) for x in fieldnames(rp)]

# Set up problem and solve it
u0 = rand(3) # initial conditions
tspan = (0.0, 300.0)
prob = ODEProblem(rossler_coupled_symbolic!, rand(6), tspan, params)
sol = solve(prob, Rosenbrock23(), saveat = 0.05)

# Plot stuff
using Plots; plotlyjs()

# 3D plot
p1 = plot(sol, vars=(1,2,3), width = 0.8, color = "black", denseplot = false)

# Plot x₁ against time (note: var 0 is time)
p2 = plot(sol, vars=(0, 1), width = 0.8, color = "black", denseplot = false)

```
"""
rossler_symbolic! = @ode_def RosslerSystem begin
    # First Rossler system
    dx = -y - z
    dy = x + a*y
    dz = b + z*(x - c)
end a b c

@with_kw struct RosslerParams
  a = 0.1 # parameter in the rossler attractor
  b = 0.1 # parameter in the rossler attractor
  c = 14.0 # parameter in the rossler attractor
end

"""
    solve_rossler(; rp = RosslerCoupledParams(),
                  tend = 300.0,
                  timestep = 0.5,
                  u0 = rand(6))

Solves the coupled rossler system (for more details, see docs for
`rossler_symbolic!``).

The solution goes from time zero until `tend` and is sampled at regular time
intervals of width `timestep`.

`u0` is a 3-element vector of initial conditions, and `rp` is an instance of
`RosslerParams`.
"""
function solve_rossler(; rp = RosslerParams(),
                        tend = 300.0,
                        timestep = 0.5,
                        u0 = rand(3))

  # Define a set of Rossler parameters. We're using the default here, but you can
  # also assign specific parameters as RosslerCoupledParams(ϵ₂ = 0.1, ω₁ = 1.1);
  # for the remaining parameters, default values are used.
  rp = RosslerCoupledParams()

  # Represent params as tuple (needed for the symbolically defined system)
  params = [getfield(rp, x) for x in fieldnames(rp)]

  # Set up problem and solve it
  tspan = (0.0, tend)
  prob = ODEProblem(rossler_symbolic!, u0, tspan, params)

  solve(prob, Rosenbrock23(), saveat = timestep)
end
