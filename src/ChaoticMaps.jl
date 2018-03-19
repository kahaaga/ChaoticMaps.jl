__precompile__(true)
module ChaoticMaps

using Reexport
########################################
# Common structures used by all the maps.
########################################
include("InformationStruct.jl")
export Info

##################################################################
# Individual maps (each map module may contain several submodules)
##################################################################
include("Logistic.jl")
export Logistic


end
