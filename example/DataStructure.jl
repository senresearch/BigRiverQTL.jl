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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# # Data Structure for genomics data

# ## Libraries

using Revise

# ### `Geno`: vector of 2D named Array of genotype data
# * one 2D array per chromosome
# * `samples` contains sample names such as genotypes or individual IDs
# * `markers` contains marker names for each chromosome
# * `val` contains alleles information

"""
`Geno` type containing genotype information for all chromosomes.

* `samples` contains sample names such as genotypes or individual IDs
* `chromosomes` contains chromosome names
* `markers` contains marker names for each chromosome
* `val` is a vector of matrices containing allele information in each chromosome
"""
struct Geno
    samples::Vector{String}
    chromosomes::Vector{String}
    markers::Vector{Vector{String}}
    val::Vector{Array{Int}}
end

?Geno

# ### `Gmap`: vector of 1D named Array of genetic map
#
# * `chromosomes` contains chromosomes names
# * `markers` contains markers's names for each chromosome
# * `val` contains relative position of markers in each chromosome

struct Gmap
    chromosomes::Vector{String}
    markers::Vector{Vector{String}}
    val::Vector{Vector{Float64}}
end


# ### pheno: vector of 2D named Array of phenotypes
# * row names are genotype
# * columns names are traits
# * value: numeric

struct Pheno
    samples::Vector{String}
    traits::Vector{String}
    val::Matrix{Float64}
end


# ### phenocov: 1D named Array of phenotype covariates
# * row names are traits
# * columns contain description
# * value: character

struct Phenocov
    traits::Vector{String}
    descriptions::Vector{String}
end


# ### pmap: vector of 1D named Array of genetic map
# * one 1D array per chromosome
# * row names are genotype
# * columns names are markers
# * value: numeric

"""
`Pmap` type contains the genetic map showing the relative location of genetic markers as phenotype.

"""
struct Pmap
    chromosomes::Vector{String}
    markers::Vector{Vector{String}}
    locations::Vector{Vector{Float64}}
    unit::String
end


# ### is_female: a named 1D-array of boolean values which individual/genotype is female

struct IsFemale
    samples::Vector{String}
    isfemale::Vector{Bool}
end


struct cross_type
    crs_typ::String
end

struct cross_info
    crs_drc:: Vector{Int64}
end

struct alleles
    allele::Vector{String}
end

struct is_X_char
    chromosomes::Vector{String}
    is_X::Vector{Bool}
end


struct BigRiverQTLData
    cross_type::cross_type
    geno::geno
    gmap::gmap
    pheno::pheno
    pmap::pmap
    phenocov::phenocov
    is_X_char::is_X_char
    is_female::is_female
    cross_info::cross_info
    alleles::alleles
end

