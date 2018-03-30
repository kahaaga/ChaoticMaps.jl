__precompile__(true)
module ChaoticMaps

using DifferentialEquations, Parameters

using Reexport

########################################
# Common structures used by all the maps.
########################################
include("InformationStruct.jl")
export
    Info

##################################################################
# Individual maps (each map module may contain several submodules)
##################################################################
include("Logistic.jl")
include("Rossler.jl")

export
    Logistic,
    Rossler
end
