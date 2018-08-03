
module TwoSpecies
using Distributions
using Parameters
using ...Info.MapInfo

"""
Parameters for the map.
"""
Parameters.@with_kw struct Params
    μXdist::Distributions.Distribution = Uniform(2, 4)
    μYdist::Distributions.Distribution = Uniform(3.7, 3.9)
    initdist::Distributions.Distribution = Uniform(0.1, 0.9)
    αdist::Distributions.Distribution = Uniform(0.02, 0.1)
    x₀::Float64 = rand(initdist)
    y₀::Float64 = rand(initdist)
    μX::Float64 = rand(μXdist)
    μY::Float64 = rand(μYdist)
    αXonY::Float64 = rand(αdist)
    αYonX::Float64 = rand(αdist)
end

"""
Logistic map object. Holds a realization of the twospecies logistic map, given by the
`itermap()` function. `X` and `Y` are time series of trajectories, `params` holds the
information used to generate the trajectiories, and `datestamp` is the time the realization
was generated.
"""
Parameters.@with_kw struct Map
    X::Vector{Float64} = Float64[]
    Y::Vector{Float64} = Float64[]
    info::MapInfo = MapInfo()
    params::Params = Params()
    datestamp::DateTime = now()
end

"""
    map(;params = LogisticMapParams)

Generate a trajectory of a bidirectionally coupled logistic map whose properties is
specified by a LogisticMapParams instance `p`.

The generated trajectory consists of two time series, X and Y, each `n_pts`
points long, starting from initial condition (x₀, y₀). If not specified otherwise, the
initial conditions are drawn according to `initdist`, which can be any
Distributions.Distribution object. The default is  `initdist= Uniform(0.01, 0.99)`.

The map is allowed to
evolve for `n_transient` time steps before we start to sample. After the
transient period is over, start sampling at every `stepsize`th iteration.

The coupling strength is given by `αXonY` and `αYonX`, which measures the
forcing strength of X on Y and the forcing strength of Y on X, in that order. If not
specified otherwise, `αXonY` and `αYonX` are drawn from `αdist`, which can be any
Distributions.Distribution object. The default is  `αdist= Uniform(0.05, 0.1)`.

The parameters `μX` and `μY` appear in the individual dynamics of X and Y, respectively.
If not specified otherwise, `μX` and `μY` are drawn from `μXdist` and `μYdist`, which both
can be any Distributions.Distribution object. The defaults are `μXdist = Uniform(2, 4)` and
`μYdist = Uniform(3.7, 3.9)`.

The `noise_frac` argument is the fraction of uniformly distributed noise relative to
1 standard deviation of the data.

If the model runs successfully given the chosen parameters, the function returns a
LogisticMap object. If the model collapses, a NaN is returned.
"""
function itermap(; p = Params(), info = MapInfo())

    X, Y = Vector{Float64}((info.n_pts*info.stepsize) + info.n_transient),
            Vector{Float64}((info.n_pts*info.stepsize) + info.n_transient)
    X[1] = p.x₀
    Y[1] = p.y₀

    # Run for `n_transient` steps, then start generating the map.
    for i = 2:info.n_transient+1
        X[i] = X[i-1] * (p.μX - p.μX*X[i-1] - p.αYonX*Y[i-1])
        Y[i] = Y[i-1] * (p.μY - p.μY*Y[i-1] - p.αXonY*X[i-1])
    end

    start_ind = info.n_transient + 1
    for i = start_ind:(info.n_pts*info.stepsize + (start_ind - 1))
        X[i] = X[i-1] * (p.μX - p.μX*X[i-1] - p.αYonX*Y[i-1])
        Y[i] = Y[i-1] * (p.μY - p.μY*Y[i-1] - p.αXonY*X[i-1])
    end


    if any(isnan, X) || any(isnan, Y)
        warn("Model collapsed. Retrying, keeping forcing strength and random seed, but randomising all other parameters.")
        return itermap(p = Params(αXonY = p.αXonY, αYonX = p.αYonX), info = info)
    end

    sample_inds = start_ind:info.stepsize:(info.n_pts*info.stepsize + (start_ind - 1))
    if info.noise_frac > 0
        X_final = X[sample_inds]
        Y_final = Y[sample_inds]
        noiseX = rand(Distributions.Uniform(-(std(X_final) * info.noise_frac), std(X_final) * info.noise_frac), info.n_pts)
        noiseY = rand(Distributions.Uniform(-(std(Y_final) * info.noise_frac), std(Y_final) * info.noise_frac), info.n_pts)
        return Map(
            X = X_final + noiseX,
            Y = Y_final + noiseY,
            params = p,
            info = info,
            datestamp = now())
    else
        return Map(X = X[sample_inds], Y = Y[sample_inds], params = p, info = info, datestamp = now())
    end
end

export Params, Map, itermap

end
