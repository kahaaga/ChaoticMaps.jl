module Info

using Parameters

@with_kw struct MapInfo
    n_pts::Int64 = 100
    n_transient::Int64 = 10000
    stepsize::Int64 = 1
    noise_frac::Float64 = 0.0
    seed::Array{Int32,1} = reinterpret(Int32, Base.Random.GLOBAL_RNG.seed)
end

export MapInfo
end
