#############################
# Generate Data for testing #
#############################

#=
gmap = Gmap(["chr1", "chr2"], [["m1", "m2"], ["m3"]], [[1.0, 2.0], [3.0]])
cross_info = CrossInfo(["sample1", "sample2"], [1, 0])
cross_type = CrossType("risib")
alleles = Alleles(["A", "B", "C"])
geno_type = GenoType(Dict("A" => 1, "H" => 2, "B" => 3))
geno_transpose = GenoTranspose(true)
geno = Geno(
	["sample1"], 
	["chr1"], 
	[["marker1"]], [[Int16(1) Int16(2); Int16(3) Int16(4)]], 
	CrossType("risib"), 
	Alleles(["A", "B"]), 
	GenoType(Dict("A" => 1)), 
	GenoTranspose(false)
)
=#
#pheno = Pheno(["individual1", "individual2"], ["trait1", "trait2"], [1.0 2.0; nothing 3.0])

########
data_dir = joinpath(@__DIR__, "data/BXD/");
file = joinpath(data_dir, "bxd.json");

# Transforming data to a optimised and accessible data type
data = get_geneticstudydata(file);

# Testing the Gmap type
gmap = data.gmap
@testset "Gmap Tests" begin

	# @test gmap.chr == ["chr1", "chr2"]
	# @test length(gmap.marker_name) == 2
	@test all(isa.(gmap.pos, Vector{Float64}))
end

# Testing the CrossType type
crosstype = data.geno.cross_type
@testset "CrossType Tests" begin
	@test isa(cross_type.type, String)
end

# Testing the CrossInfo type
cross_info = data.cross_info
@testset "CrossInfo Tests" begin
	@test length(cross_info.sample_id) == length(cross_info.direction)
	@test all(isa.(cross_info.direction, Int))
end

# Testing the Alleles type
alleles = data.geno.alleles
@testset "Alleles Tests" begin
	# @test length(alleles.val) == 3
	@test all(isa.(alleles.val, String))
end

# Testing the GenoType type
#@testset "GenoType Tests" begin
# @test geno_type.label["A"] == 1
# @test geno_type.label["H"] == 2
#end

# Testing the GenoTranspose type
#@testset "GenoTranspose Tests" begin
# @test geno_transpose.val == true
#end

# Testing the Geno type
#@testset "Geno Tests" begin
# @test length(geno.sample_id) == 1
#  @test size(geno.val[1]) == (2, 2)
#end

# Testing the Pheno type
pheno = data.pheno
@testset "Pheno Tests" begin
	@test size(pheno.val, 1) == length(pheno.sample_id)
	# @test pheno.val[2, 1] === nothing
end


# Checks whether the output of `get_geneticstudydata` function is of type `GeneticStudyData`
function BigRiverQTLData_struct_test(filename::String, testname::String)
	@testset "BigRiverQTLData_struct_test" begin

		println("Test results for checking whether the output of `get_geneticstudydata` function is of type `GeneticStudyData`: ",
			@test typeof(get_geneticstudydata(filename)) == GeneticStudyData)


	end

end




#########################
# Test: BigRiverQTLData #
#########################

data_dir = joinpath(@__DIR__, "data/BXD/");
file = joinpath(data_dir, "bxd.json");


BigRiverQTLData_struct_test(file, "get_bigriverqtldata")

