#############################
# Generate Data for testing #
#############################

data_dir = joinpath(@__DIR__, "../data/BXD/")
file = joinpath(data_dir, "bxd.json")



# Checks whether the output of `get_bigriverqtldata` function is of type `BigRiverQTLData`
function BigRiverQTLData_struct_test(filename::String, testname::String)
    @testset "$testname" begin

        println("Test results for checking whether the output of `get_bigriverqtldata` function is of type `BigRiverQTLData`: ", @test typeof(get_bigriverqtldata(filename))==BigRiverQTLData)
         
        
    end
        
end




#########################
# Test: BigRiverQTLData #
#########################



BigRiverQTLData_struct_test(file, "get_bigriverqtldata")

