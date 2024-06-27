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

# ### geno: vector of 2D named Array of genotype data
# * one 2D array per chromosome
# * row names are genotype
# * columns names are markers
# * value: numeric

struct geno2
    gno_name::Vector{String}
    chr_name::Vector{String}
    mkrs_name::Vector{Vector{String}}
    val::Vector{Array{Float64}}
end

# ### gmap: vector of 1D named Array of genetic map
# * one 1D array per chromosome
# * row names are genotype
# * columns names are markers
# * value: numeric

struct gmap
    chr_name::Vector{String}
    mkrs_name::Vector{Vector{String}}
    val::Vector{Vector{Float64}}
end
    

# ### pheno: vector of 2D named Array of phenotypes
# * row names are genotype
# * columns names are traits
# * value: numeric

struct pheno
    gno_name::Vector{String}
    traits::Vector{String}
    val::Array{Float64,2}
end
    

# ### phenocov: 1D named Array of phenotype covariates
# * row names are traits
# * columns contain description
# * value: character

struct phenocov
    traits::Vector{String}
    dscrpt::Vector{String}
end
    

# ### pmap: vector of 1D named Array of genetic map
# * one 1D array per chromosome
# * row names are genotype
# * columns names are markers
# * value: numeric

struct pmap
    chr_name::Vector{String}
    mkrs_name::Vector{Vector{String}}
    val::Vector{Vector{Float64}}
end
    

# ### is_female: a named 1D-array of boolean values which individual/genotype is female

struct is_female
    chr_name::Vector{String}
    is_f::Vector{Bool}
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
    is_X::Vector{Bool}
end

 struct geno_data1
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
    
