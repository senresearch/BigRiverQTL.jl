
"""
`Gmap` type contains markers names and their relative positions in each chromosome.

* `chr` contains chromosomes names.
* `marker_name` contains marker name's names for each chromosome.
* `pos` is a vector of vector  containing relative position of marker_name in each chromosome.
"""
struct Gmap
	chr::Vector{AbstractString}
	marker_name::Vector{Vector{AbstractString}}
	pos::Vector{Vector{Float64}}
end


"""
`CrossType` type contains the cross type, for example: risib => Recombinant inbred lines (RILs) by sibling mating.
* `type` is a string indicating the type of the cross.
"""
struct CrossType
	type::AbstractString
end


"""
`CrossInfo` type contains information about the cross direction of sample_id.
* `sample_id` contains sample names such as genotypes or individual IDs.
* `direction` is a vector containing the cross direction of sample_id.
"""
struct CrossInfo
	sample_id::Vector{AbstractString}
	direction::Vector{AbstractString}
end


"""
`Alleles` type contains the names of the alleles.
* `val` is a vector containing the names of the alleles.
"""
struct Alleles
	val::Vector{AbstractString}
end


"""
`GenoType` type indicates the transformation or label for genotypes. For example, "(A=1, H=2, B=3, D=4, C=5)"
* `label` is a dictionary type indicating the transformation or label for genotypes.
"""
struct GenoType
	label::Dict
end


"""
`GenoTranspose` type indicates whether the geno matrix is transposed.
* `val` is a boolean type indicating whether the geno matrix is transposed.
"""
struct GenoTranspose
	val::Bool
end


"""
`Geno` type contains genotype information for all chromosomes.

* `sample_id` contains sample names such as genotypes or individual IDs.
* `chromosomes` contains chromosome names.
* `marker_name` contains marker names for each chromosome.
* `val` is a vector of matrices containing allele information in each chromosome.
"""
struct Geno
	sample_id::Vector{AbstractString}
	chromosomes::Vector{AbstractString}
	marker_name::Vector{Vector{AbstractString}}
	val::Vector{Matrix{Union{Missing, Int16}}}
	cross_type::CrossType
	alleles::Alleles
	geno_type::GenoType
	geno_transpose::GenoTranspose

end


"""
 `Pmap` type contains the genetic map showing the relative location of genetic marker_name as phenotype.
* `chromosomes` contains chromosomes names.
* `marker_name` contains marker  names for each chromosome.
* `pos` is a vector of vector containing relative position of marker_name as phenotypes in each chromosome.
* `unit` contains unit for the chromosome length.
"""
struct Pmap
	chromosomes::Vector{AbstractString}
	marker_name::Vector{Vector{AbstractString}}
	pos::Vector{Vector{Float64}}
	unit::AbstractString
end


"""
`Pheno` type contains phenotypes data.
* `sample_id` contains sample names such as genotypes or individual IDs.
* `traits` contains trait names.
*  `val` is a matrix containing phenotype/ traits values.
"""
struct Pheno
	sample_id::Vector{AbstractString}
	traits::Vector{AbstractString}
	val::Matrix{Union{Nothing, Float64}}
end


"""
`Phenocov` type contains the description of the phenotypes.
* `traits` contains trait names.
* `descriptions` is a vector containing the description for each phenotype.
"""
struct Phenocov
	traits::Vector{AbstractString}
	descriptions::Vector{AbstractString}
end


"""
IsFemale type indicates if the sample_id (genotypes or individuals) are females.
* `sample_id` contains sample names such as genotypes or individual IDs.
* `val` is a vector containing boolean values indicating if each sample (genotype or individual) is a female.
"""
struct IsFemale
	sample_id::Vector{AbstractString}
	val::Vector{Bool}
end


"""
`IsXChar` type indicates which chromosome is the X one.
* `chromosomes` contains chromosome names.
* `val` is a vector of boolean values indicating which chromosome is the X one.
"""
struct IsXChar
	chromosomes::Vector{AbstractString}
	val::Vector{Bool}
end


"""
`GeneticStudyData` type contains genomics stuctured data suitable to use for QTL analysis.
* `geno` is a field of type `Geno`. Refer to `Geno` type for more imformation.
* `gmap` is a field of type `Gmap`. Refer to `Gmap` type for more imformation.
* `pheno` is a field of type `Pheno`. Refer to `Pheno` type for more imformation.
* `pmap` is a field of type `Pmap`. Refer to `Pmap` type for more imformation.
* `phenocov` is a field of type `Phenocov`. Refer to `Phenocov` type for more imformation.
* `isXchar` is a field of type `IsXChar`. Refer to `IsXChar` type for more imformation.
* `isfemale` is a field of type `IsFemale`. Refer to `IsFemale` type for more imformation.
* `crosstype` is a field of type `CrossType`. Refer to `CrossType` type for more imformation.
* `crossinfo` is a field of type `CrossInfo`. Refer to `CrossInfo` type for more imformation.
* `alleles` is a field of type `Alleles`. Refer to `Alleles` type for more imformation.
* `genotype` is a field of type `Genotype`. Refer to `Genotype` type for more imformation.
* `genotranspose` is a field of type `GenoTranspose`. Refer to `GenoTranspose` type for more imformation.
"""
struct GeneticStudyData
	gmap::Gmap
	geno::Geno
	pmap::Pmap
	pheno::Pheno
	phenocov::Phenocov
	isXchar::IsXChar
	isfemale::IsFemale
	crossinfo::CrossInfo
end
