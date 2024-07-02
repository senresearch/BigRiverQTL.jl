
#############################
# Generate Data for testing #
#############################

data_file = joinpath(@__DIR__, "data/geno_proba.csv")
geno1 = readdlm(data_file, ',', Float64, '\n')

# Another constructed data set
geno2 =  [1    2    3    4    5    6; 
          1    0    0    0    0    0; 
          1001 1002 1003 1004 1005 1006; 
          1    2    3    4    5    6.1] |> x -> x ./ sum(x, dims= 2);



#=
* Dimensions: $n \times n$ where $n$ is the number of subjects.
* All elements in the diagonals are $1$.
* All elements in each kinship matrix lie between $-1$ and $1$.
* Symmetric
* Positive Definite
=#







# Checks the above mentioned properties
function kinship_test(kinmat::Matrix{Float64}, testname::String)
    @testset "$testname" begin

        println("Test kinship dimensions: ", @test size(kinmat,1)==size(kinmat,2))
        println("Test range: ", @test (maximum(kinmat)<=1.0) && (minimum(kinmat)>=-1.0))
        println("Test diagonal: ", @test diag(kinmat) ≈ ones(size(kinmat,1)))  
        println("Test diagonal: ", @test isposdef(kinmat))  
        
    end
        
end

#####################
# Test: calckinship #
#####################

K1 = calckinship(geno1)
K2 = calckinship(geno2)

kinship_test(K1, "calckinship()")
