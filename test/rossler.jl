using ChaoticMaps
using ChaoticMaps.Rossler

@testset "Rossler Systems" begin
    tend = 100.0
    timestep = 1.0

    @testset "Single Rossler attractor" begin
        u0 = rand(3) # initial conditions
        params = RosslerParams()
        sol = solve_rossler(rp = params, tend = tend, timestep = timestep, u0 = u0)
        @test length(sol) == 101
    end

    @testset "Two coupled Rossler attractors" begin
        u0 = rand(6) # initial conditions
        p1 = RosslerCoupledParams(ϵ₁ = 0.1, ϵ₂ = 0.0) # 1 forces 2
        p2 = RosslerCoupledParams(ϵ₁ = 0.0, ϵ₂ = 0.1) # 2 forces 1
        p3 = RosslerCoupledParams(ϵ₁ = 0, ϵ₂ = 0)     # no forcing

        sol1 = solve_rossler_coupled(rp = p1, tend = tend, timestep = timestep, u0 = u0)
        sol2 = solve_rossler_coupled(rp = p2, tend = tend, timestep = timestep, u0 = u0)
        sol3 = solve_rossler_coupled(rp = p3, tend = tend, timestep = timestep, u0 = u0)

        @test length(sol1.t) == 101
        @test length(sol2.t) == 101
        @test length(sol3.t) == 101
    end

end
