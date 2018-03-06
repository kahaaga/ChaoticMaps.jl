module Logistic

using Distributions

""" Logistic map parameters. """
type LogisticMapParams
    μXdist::Distributions.Distribution
    μYdist::Distributions.Distribution
    initdist::Distributions.Distribution
    αdist::Distributions.Distribution
    x₀::Float64
    y₀::Float64
    μX::Float64
    μY::Float64
    αXonY::Float64
    αYonX::Float64
    n_pts::Int
    n_transient::Int
    stepsize::Int
    seed::Any
end

""" Logistic map object. Holds a realization of the logistic map, given by the
`logistic_map()` function. `X` and `Y` are time series of trajectories, `params` holds the
information used to generate the trajectiories, and `datestamp` is the time the realization
was generated."""
type LogisticMap
    X::Vector{Float64}
    Y::Vector{Float64}
    params::LogisticMapParams
    datestamp::DateTime
end

"""
    logistic_map(;
        μXdist::Distributions.Distribution = Uniform(2, 4),
        μYdist::Distributions.Distribution = Uniform(3.7, 3.9),
        initdist::Distributions.Distribution = Uniform(0.01, 0.99),
        αdist::Distributions.Distribution = Uniform(0.05, 0.1),
        x₀::Float64 = rand(initdist),
        y₀::Float64 = rand(initdist),
        μX::Float64 = rand(μXdist),
        μY::Float64 = rand(μYdist),
        αXonY::Float64 = rand(αdist),
        αYonX::Float64 = rand(αdist),
        n_pts::Int = 100,
        n_transient::Int = 10000,
        stepsize::Int = 1,
        noise_frac::Float64 = 0.0,
        seed::Int = 1234)

Generate a trajectory of a bidirectionally coupled logistic map.

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
function logistic_map(;
    μXdist::Distributions.Distribution = Uniform(2, 4),
    μYdist::Distributions.Distribution = Uniform(3.7, 3.9),
    initdist::Distributions.Distribution = Uniform(0.01, 0.99),
    αdist::Distributions.Distribution = Uniform(0.05, 0.1),
    x₀::Float64 = rand(initdist),
    y₀::Float64 = rand(initdist),
    μX::Float64 = rand(μXdist),
    μY::Float64 = rand(μYdist),
    αXonY::Float64 = rand(αdist),
    αYonX::Float64 = rand(αdist),
    n_pts::Int = 100,
    n_transient::Int = 10000,
    stepsize::Int = 1,
    noise_frac::Float64 = 0.0,
    seed::Int = 1234)

    X, Y = Vector{Float64}((n_pts*stepsize) + n_transient),
            Vector{Float64}((n_pts*stepsize) + n_transient)
    X[1] = x₀
    Y[1] = y₀

    # Run for `n_transient` steps, then start generating the map.
    for i = 2:n_transient+1
        X[i] = X[i-1] * (μX - μX*X[i-1] - αYonX*Y[i-1])
        Y[i] = Y[i-1] * (μY - μY*Y[i-1] - αXonY*X[i-1])
    end

    start_ind = n_transient + 1
    for i = start_ind:(n_pts*stepsize + (start_ind - 1))
        X[i] = X[i-1] * (μX - μX*X[i-1] - αYonX*Y[i-1])
        Y[i] = Y[i-1] * (μY - μY*Y[i-1] - αXonY*X[i-1])
    end

	params = LogisticMapParams(
                μXdist, μYdist, initdist, αdist,
                x₀, y₀,
                μX, μY,
                αXonY, αYonX,
                n_pts, n_transient, stepsize,
                reinterpret(Int32, Base.Random.GLOBAL_RNG.seed)
            )

     if any(isnan, X) || any(isnan, Y)
        warn("Model collapsed. Returning NaN")
        return NaN
     end

    sample_inds = start_ind:stepsize:(n_pts*stepsize + (start_ind - 1))
    if noise_frac > 0
        X_final = X[sample_inds]
        Y_final = Y[sample_inds]
        noiseX = rand(Uniform(-(std(X_final) * noise_frac), std(X_final) * noise_frac), n_pts)
        noiseY = rand(Uniform(-(std(Y_final) * noise_frac), std(Y_final) * noise_frac), n_pts)
        return LogisticMap(X_final + noiseX, Y_final + noiseY, params, now())
    else
        return LogisticMap(X[sample_inds], Y[sample_inds], params, now())
    end
end


"""
Exactly the same as the function `logistic_map`, but fixing the forcing strengths
between the components.
"""
function logistic_sample_fixed_couplingstrength(
    αXonY::Float64,
    αYonX::Float64;
    n_pts::Int = 100,
    n_transient::Int = 10000,
    stepsize::Int = 1,
    noise_frac::Float64 = 0.0,
    seed::Int = 10)

    logistic_map(
        n_pts = n_pts,
        n_transient = n_transient,
        stepsize = stepsize,
        noise_frac = noise_frac,
        seed = seed
    )
end

# Trigger precompilation.
l = logistic_sample_fixed_couplingstrength(0.1, 0.3, n_pts = 1000, noise_frac = 0.05)
end
