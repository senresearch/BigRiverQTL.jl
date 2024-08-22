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


#########################
# Testing the Gmap type #
#########################
@testset "Gmap Tests" begin
	gmap = get_gmap(file)
	@test all(isa.(gmap.pos, Vector{Float64}))
end


##############################
# Testing the CrossType type #
##############################
@testset "CrossType Tests" begin
	cross_type = get_crosstype(file)
	@test isa(cross_type.type, String)
end


############################
# Testing the Alleles type #
############################
@testset "Alleles Tests" begin
	alleles = get_alleles(file)
	# @test length(alleles.val) == 3
	@test all(isa.(alleles.val, String))
end


#############################
# Testing the GenoType type #
#############################
@testset "GenoType Tests" begin
	genotype_dict = get_genotype(file)
	@test genotype_dict.label["B"] == 1
	@test genotype_dict.label["D"] == 2
end


##################################
# Testing the GenoTranspose type #
##################################
@testset "GenoTranspose Tests" begin
	geno_transpose = get_genotranspose(file)
	@test geno_transpose.val == true
end


#########################
# Testing the Geno type #
#########################
@testset "Geno Tests" begin
	geno = get_geno(file)
	@test length(geno.sample_id) == 198
	@test size(geno.val[1]) == (198, 636)
end


#########################
# Testing the Pmap type #
#########################
@testset "Pmap Tests" begin
	pmap = get_pmap(file)
	@test length(pmap.pos[1]) == 636
	@test pmap.unit == "Mbp"
end


#########################
# Testing the Phenotype #
#########################
@testset "Pheno Tests" begin
	pheno = get_pheno(file)
	@test length(pheno.sample_id) == 198
	@test length(pheno.traits) == 5806
	@test isapprox(pheno.val[1:2, 1:2], [61.4  54.1;49.0  50.1], atol=1e-4)
end


#############################
# Testing the Phenocov type #
#############################
@testset "Phenocov Tests" begin
	phenocovar = get_phenocovar(file)
	@test length(phenocovar.traits) == 5806
	@test phenocovar.descriptions[5800] == "MFL-54"
end


##############################
# Testing the CrossInfo type #
##############################
@testset "CrossInfo Tests" begin
	cross_info = get_crossinfo(file)
	@test length(cross_info.sample_id) == 198
	@test cross_info.direction[1] == "BxD"
end


############################
# Testing the IsXChar type #
############################
@testset "IsXChar Tests" begin
	isxchar = get_isxchar(file)
	@test length(isxchar.chr) == 20
	@test isxchar.val[20] == true
end


#############################
# Testing the IsFemale type #
#############################
@testset "IsFemale Tests" begin
	isfemale = get_isfemale(file)
	@test length(isfemale.sample_id) == 198
	@test isfemale.val[1] == true
end


#####################################
# Testing the GeneticStudyData type #
#####################################
@testset "GeneticStudyData Test" begin
	geneticstudydata = get_geneticstudydata(file)
	@test isa(geneticstudydata, GeneticStudyData)
end

