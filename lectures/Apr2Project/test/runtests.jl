using Test
using Apr2Project

@testset "Test add function" begin
    @test include("functiontest1.jl")
    @test include("functiontest2.jl")
end

@testset "Test subtract function" begin
    @test include("subtest1.jl")
end
