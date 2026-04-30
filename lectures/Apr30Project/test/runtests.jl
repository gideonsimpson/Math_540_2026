using Test
using Apr30Project
using LinearAlgebra

@testset "Sparse Differentiation" begin
    @test include("sparse1.jl")
end
