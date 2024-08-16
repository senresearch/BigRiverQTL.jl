###########
# Library #
###########

using BigRiverQTL, LinearAlgebra
using Statistics
using Test
using DelimitedFiles, Plots

########
# Test #
########

@testset "BigRiverQTL" begin 
    include("kinship_test.jl")
    include("datastruc_test.jl")
    include("utils_test.jl")
    include("plots_test.jl") 
end