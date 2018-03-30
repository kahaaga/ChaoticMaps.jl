module Rossler

using Parameters, DifferentialEquations

include("rossler/rossler.jl")
include("rossler/rossler_coupled.jl")

export
    RosslerParams, rossler_symbolic!,solve_rossler,
    RosslerCoupledParams, rossler_coupled_symbolic!, solve_rossler_coupled
end
