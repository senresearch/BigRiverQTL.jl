using Revise

"""
`Geno` type

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

struct Gmap
    chromosomes::Vector{String}
    markers::Vector{Vector{String}}
    val::Vector{Vector{Float64}}
end

struct pheno
    gno_name::Vector{String}
    traits::Vector{String}
    val::Array{Float64,2}
end

struct phenocov
    traits::Vector{String}
    dscrpt::Vector{String}
end

struct pmap
    chr_name::Vector{String}
    mkrs_name::Vector{Vector{String}}
    val::Vector{Vector{Float64}}
end

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
