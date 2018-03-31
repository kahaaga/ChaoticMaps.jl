using DifferentialEquations, Parameters

"""
Parameters for the Lorenz system.
"""
@with_kw struct LorenzParams
    σ::Float64 = 10
    ρ::Float64 = 28
    β::Float64 = 8/3
end


"""
    lorenz!(du, u, p, t)

Lorentz system.

# Examples

```julia-repl
# Set up an ODEProblem for the Lorenz system.
u0 = [1.0, 0.0, 0.0]
params = LorenzParams()
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan, params)

# Solve the system.
sol = solve(prob)

# Plot the solution
using Plots; plotlyjs()
plot(sol, vars=(1,2,3), width = 1.2, opacity = 0.8, color = "black")
```
"""
function lorenz!(du, u, p, t)
    du[1] = p[1]*(u[2] - u[1])
    du[2] = u[1]*(p[2] - u[3]) - u[2]
    du[3] = u[1]*u[2] - p[3]*u[3]
end

"""
  lorenz!(du, u, p::LorentzParams, t)

Lorentz system with parameter struct.

# Examples

```julia-repl
# Set up an ODEProblem for the Lorenz system.
u0 = [1.0, 0.0, 0.0]
params = LorenzParams(σ = 10, ρ = 28, β = 8/3)
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan, params)

# Solve the system.
sol = solve(prob)

# Plot the solution
using Plots; plotlyjs()
plot(sol, vars=(1,2,3), width = 1.2, opacity = 0.8, color = "black")
```
"""
function lorenz!(du, u, p::LorenzParams, t)
    du[1] = p.σ*(u[2] - u[1])
    du[2] = u[1]*(p.ρ - u[3]) - u[2]
    du[3] = u[1]*u[2] - p.β*u[3]
end


"""
  lorenz_symbolic! = @ode_def Lorenz begin
      dx = σ*(x - y)
      dy = x*(ρ - z) - y
      dz = x*y - β*z
  end σ ρ β

Lorentz system from symbolic expression. Parameters must be fed as a tuple.
See example.

# Examples

```julia-repl
# Set up an ODEProblem for the Lorenz system.
u0 = [1.0, 0.0, 0.0]
params = (10, 28, 8/3) # (σ, ρ, β)
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz_symbolic!, rand(3), tspan, params)

# Solve the system.
sol = solve(prob)

# Plot the solution
using Plots; plotlyjs()
plot(sol, vars=(1,2,3), width = 1.2, opacity = 0.8, color = "black")

"""
lorenz_symbolic! = @ode_def LorenzSystem begin
    dx = σ*(y - x)
    dy = x*(ρ - z) - y
    dz = x*y - β*z
end σ ρ β


"""
    solve_lorenz(; p = LorenzParams(),
                  tend = 300.0,
                  timestep = 0.5,
                  u0 = rand(3))

Solves the coupled rossler system (for more details, see docs for
`lorenz_symbolic!`).

The solution goes from time zero until `tend` and is sampled at regular time
intervals of width `timestep`.

`u0` is a 3-element vector of initial conditions, and `p` is an instance of
`LorenzParams`.
"""
function solve_lorenz(; p = LorenzParams(),
                        tend = 300.0,
                        timestep = 0.5,
                        u0 = rand(3))

  # Define a set of Lorenz parameters. We're using the default here.
  # Represent params as tuple (needed for the symbolically defined system)
  params = [getfield(p, x) for x in fieldnames(p)]

  # Set up problem and solve it
  tspan = (0.0, tend)

  prob = ODEProblem(lorenz_symbolic!, u0, tspan, params)

  solve(prob, Rosenbrock23(), saveat = timestep)
end
