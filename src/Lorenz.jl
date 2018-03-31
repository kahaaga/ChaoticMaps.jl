module Lorenz


using Parameters, DifferentialEquations

include("lorenz/lorenz.jl")

export
    LorenzParams, lorenz!, lorenz_symbolic!, solve_lorenz
end
