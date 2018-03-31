using ChaoticMaps
using ChaoticMaps.Lorenz

@testset "Lorenz Systems" begin
    tend = 100.0
    timestep = 1.0

    @testset "Single Lorenz attractor" begin
        u0 = rand(3) # initial conditions
        params = LorenzParams()
        sol = solve_lorenz(p = params, tend = tend, timestep = timestep, u0 = u0)
        @test length(sol) == 101
    end
end
