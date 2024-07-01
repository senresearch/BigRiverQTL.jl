###########
# Library #
###########

using BigRiverQTL, LinearAlgebra,
using DataFrames 
using Random, Statistics
using Test


########
# Test #
########


@testset "BigRiverQTL" begin 
    include("kinship_test.jl")
    include("datastruct_test.jl")
    
end