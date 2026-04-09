using Test
using Apr9Project

@testset "Explicit Euler" begin
    @test include("fwd1.jl")
    @test include("fwd2.jl")
    @test include("fwd3.jl")
    @test include("fwd4.jl")
end

@testset "Implicit Euler" begin
    @test include("bwd1.jl")
    @test include("bwd2.jl")
end
