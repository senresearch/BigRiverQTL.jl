# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Julia_6_Threads 1.10.3
#     language: julia
#     name: julia_6_threads-1.10
# ---

# # Data Structure for genomics data

# ## Libraries

using Revise

# ### `Geno` type containing genotype information for all chromosomes.
# * `samples` contains sample names such as genotypes or individual IDs.
# * `chromosomes` contains chromosome names.
# * `markers` contains marker names for each chromosome.
# * `val` is a vector of matrices containing allele information in each chromosome.

"""
`Geno` type containing genotype information for all chromosomes.

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

# ### `Gmap` type containing relative position of markers in each chromosome.
#
# * `chromosomes` contains chromosomes names.
# * `markers` contains markers's names for each chromosome.
# * `val` is a vector of vector  containing relative position of markers in each chromosome.

"""
`Gmap` type containing relative position of markers in each chromosome.

* `chromosomes` contains chromosomes names.
* `markers` contains markers's names for each chromosome.
* `val` is a vector of vector  containing relative position of markers in each chromosome.
"""
struct Gmap
    chromosomes::Vector{String}
    markers::Vector{Vector{String}}
    val::Vector{Vector{Float64}}
end


# ### `Pheno` type containing phenotypes information.
# * `samples` contains sample names such as genotypes or individual IDs.
# * `traits` contains trait names.
# *  `val` is a matrix containing phenotype/ traits values.  

"""
`Pheno` type containing phenotypes information.
* `samples` contains sample names such as genotypes or individual IDs.
* `traits` contains trait names.
*  `val` is a matrix containing phenotype/ traits values.
"""
struct Pheno
    samples::Vector{String}
    traits::Vector{String}
    val::Matrix{Float64}
end


# ### `Phenocov`: type containing the description of the phenotypes.
# * `traits` contains trait names.
# * `descriptions` is a vector containing the description for each phenotype.
#   

"""
`Phenocov`: type containing the description of the phenotypes.
* `traits` contains trait names.
* `descriptions` is a vector containing the description for each phenotype.
"""
struct Phenocov
    traits::Vector{String}
    descriptions::Vector{String}
end


# ### `Pmap`: type contains the genetic map showing the relative location of genetic markers as phenotype.
# * `chromosomes` contains chromosomes names
# * `markers` contains markers's names for each chromosome
# * `val` is a vector of vector containing relative position of markers as phenotypes in each chromosome

"""
 `Pmap`: type contains the genetic map showing the relative location of genetic markers as phenotype.
* `chromosomes` contains chromosomes names.
* `markers` contains markers's names for each chromosome.
* `val` is a vector of vector containing relative position of markers as phenotypes in each chromosome.
"""
struct Pmap
    chromosomes::Vector{String}
    markers::Vector{Vector{String}}
    locations::Vector{Vector{Float64}}
    unit::String
end


# ### IsFemale: type indicating if the samples (genotypes or individuals) are females.
# * `samples` contains sample names such as genotypes or individual IDs.
# * `isfemale` is a vector containing boolean values indicating if each sample (genotype or individual) is a female.

"""
IsFemale: type indicating if the samples (genotypes or individuals) are females.
* `samples` contains sample names such as genotypes or individual IDs.
* `isfemale` is a vector containing boolean values indicating if each sample (genotype or individual) is a female.
"""
struct IsFemale
    samples::Vector{String}
    isfemale::Vector{Bool}
end


# ### `CrossType`: type containing the cross type for example: risib => Recombinant inbred lines (RILs) by sibling mating.
# * crs_typ is a string indicating the type of the cross.

"""
`CrossType`: type containing the cross type for example: risib => Recombinant inbred lines (RILs) by sibling mating.
* crs_typ is a string indicating the type of the cross.
"""
struct CrossType
    crs_typ::String
end

# ### `CrossInfo`: type containing information about the cross direction of samples.
# * `samples` contains sample names such as genotypes or individual IDs.
# * `crsdrc` is a vector containing the cross direction of samples.

"""
`CrossInfo`: type containing information about the cross direction of samples.
* `samples` contains sample names such as genotypes or individual IDs.
* `crsdrc` is a vector containing the cross direction of samples.
"""
struct CrossInfo
    samples::Vector{String}
    crsdrc:: Vector{Int}
end

# ### `Alleles`: type containing the names of the alleles.
# * `allele` is a vector containing the names of the alleles.

"""
`Alleles`: type containing the names of the alleles.
* `allele` is a vector containing the names of the alleles.
"""
struct Alleles
    allele::Vector{String}
end

# ### `IsXChar`: type indicating which chromosome is the X one.
# * `chromosomes` contains chromosome names.
# * `isX` is a vector of boolean values indicating which chromosome is the X one.

"""
`IsXChar`: type indicating which chromosome is the X one.
* `chromosomes` contains chromosome names.
* `isX` is a vector of boolean values indicating which chromosome is the X one.
"""
struct IsXChar
    chromosomes::Vector{String}
    isX::Vector{Bool}
end


# ### `BigRiverQTLData`: type containing genomics data suitable to use for QTL analysis using the BigRiverQTl package.
# * `geno` is a object of type `Geno`. Refer to `Geno` type for more imformation.
# * `gmap` is a object of type `Gmap`. Refer to `Gmap` type for more imformation.
# * `pheno` is a object of type `Pheno`. Refer to `Pheno` type for more imformation.
# * `pmap` is a object of type `Pmap`. Refer to `Pmap` type for more imformation.
# * `phenocov` is a object of type `Phenocov`. Refer to `Phenocov` type for more imformation.
# * `isXchar` is a object of type `IsXChar`. Refer to `IsXChar` type for more imformation.
# * `isfemale` is a object of type `IsFemale`. Refer to `IsFemale` type for more imformation.
# * `crosstype` is a object of type `CrossType`. Refer to `CrossType` type for more imformation.
# * `crossinfo` is a object of type `CrossInfo`. Refer to `CrossInfo` type for more imformation.
# * `alleles` is a object of type `Alleles`. Refer to `Alleles` type for more imformation.

"""
`BigRiverQTLData`: type containing genomics data suitable to use for QTL analysis using the BigRiverQTl package.
* `geno` is a object of type `Geno`. Refer to `Geno` type for more imformation.
* `gmap` is a object of type `Gmap`. Refer to `Gmap` type for more imformation.
* `pheno` is a object of type `Pheno`. Refer to `Pheno` type for more imformation.
* `pmap` is a object of type `Pmap`. Refer to `Pmap` type for more imformation.
* `phenocov` is a object of type `Phenocov`. Refer to `Phenocov` type for more imformation.
* `isXchar` is a object of type `IsXChar`. Refer to `IsXChar` type for more imformation.
* `isfemale` is a object of type `IsFemale`. Refer to `IsFemale` type for more imformation.
* `crosstype` is a object of type `CrossType`. Refer to `CrossType` type for more imformation.
* `crossinfo` is a object of type `CrossInfo`. Refer to `CrossInfo` type for more imformation.
* `alleles` is a object of type `Alleles`. Refer to `Alleles` type for more imformation.
"""
struct BigRiverQTLData
    geno::Geno
    gmap::Gmap
    pheno::Pheno
    pmap::Pmap
    phenocov::Phenocov
    isXchar::IsXChar
    isfemale::IsFemale
    crosstype::CrossType
    crossinfo::CrossInfo
    alleles::Alleles
end
