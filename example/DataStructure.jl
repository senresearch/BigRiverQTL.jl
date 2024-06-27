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

struct geno2
    gno_name::Vector{String}
    chr_name::Vector{String}
    mkrs_name::Vector{Vector{String}}
    val::Vector{Array{Float64}}
end

struct is_female
    is_f::Vector{Bool}
end


IS_F=is_female([true, false, true])

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

struct pheno
    gno_name::Vector{String}
    traits::Vector{String}
    val::Array{Float64,2}
end
    

struct phenocov
    traits::Vector{String}
    dscrpt::Vector{String}
end
    

struct gmap
    chr_name::Vector{String}
    mkrs_name::Vector{Vector{String}}
    val::Vector{Vector{Float64}}
end
    

struct pmap
    chr_name::Vector{String}
    mkrs_name::Vector{Vector{String}}
    val::Vector{Vector{Float64}}
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
    
