
"""
`Gmap` type contains relative position of markers in each chromosome.

* `chr` contains chromosomes names.
* `marker` contains markers's names for each chromosome.
* `pos` is a vector of vector  containing relative position of markers in each chromosome.
"""
struct Gmap
    chr::Vector{String}
    marker::Vector{Vector{String}}
    pos::Vector{Vector{Float64}}
end

"""
`Geno` type contains genotype information for all chromosomes.

* `samples` contains sample names such as genotypes or individual IDs.
* `chromosomes` contains chromosome names.
* `markers` contains marker names for each chromosome.
* `val` is a vector of matrices containing allele information in each chromosome.
"""
struct Geno
    samples::Vector{String}
    chromosomes::Vector{String}
    markers::Vector{Vector{String}}
    val::Vector{Array{Int}}
end

"""
 `Pmap` type contains the genetic map showing the relative location of genetic markers as phenotype.
* `chromosomes` contains chromosomes names.
* `markers` contains markers's names for each chromosome.
* `pos` is a vector of vector containing relative position of markers as phenotypes in each chromosome.
* `unit` contains unit for the chromosome length.
"""
struct Pmap
    chromosomes::Vector{String}
    markers::Vector{Vector{String}}
    pos::Vector{Vector{Float64}}
    unit::String
end

"""
`Pheno` type contains phenotypes data.
* `samples` contains sample names such as genotypes or individual IDs.
* `traits` contains trait names.
*  `val` is a matrix containing phenotype/ traits values.
"""
struct Pheno
    samples::Vector{String}
    traits::Vector{String}
    val::Matrix{Float64}
end

"""
`Phenocov` type contains the description of the phenotypes.
* `traits` contains trait names.
* `descriptions` is a vector containing the description for each phenotype.
"""
struct Phenocov
    traits::Vector{String}
    descriptions::Vector{String}
end


"""
IsFemale type indicates if the samples (genotypes or individuals) are females.
* `samples` contains sample names such as genotypes or individual IDs.
* `val` is a vector containing boolean values indicating if each sample (genotype or individual) is a female.
"""
struct IsFemale
    samples::Vector{String}
    val::Vector{Bool}
end

"""
`CrossType` type contains the cross type, for example: risib => Recombinant inbred lines (RILs) by sibling mating.
* `type` is a string indicating the type of the cross.
"""
struct CrossType
    type::String
end

"""
`CrossInfo` type contains information about the cross direction of samples.
* `samples` contains sample names such as genotypes or individual IDs.
* `direction` is a vector containing the cross direction of samples.
"""
struct CrossInfo
    samples::Vector{String}
    direction:: Vector{Int}
end

"""
`Alleles` type contains the names of the alleles.
* `val` is a vector containing the names of the alleles.
"""
struct Alleles
    val::Vector{String}
end

"""
`IsXChar` type indicates which chromosome is the X one.
* `chromosomes` contains chromosome names.
* `val` is a vector of boolean values indicating which chromosome is the X one.
"""
struct IsXChar
    chromosomes::Vector{String}
    val::Vector{Bool}
end

"""
`BigRiverQTLData` type contains genomics stuctured data suitable to use for QTL analysis.
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
"""
struct BigRiverQTLData
    gmap::Gmap
    geno::Geno
    pmap::Pmap
    pheno::Pheno
    phenocov::Phenocov
    isXchar::IsXChar
    isfemale::IsFemale
    crosstype::CrossType
    crossinfo::CrossInfo
    alleles::Alleles
end
