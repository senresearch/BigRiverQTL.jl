"""
loco_processing(geno::DataFrame, genome_start_index::Int)

Takes in a genotype DataFrame and the index of the first actual genome data column (after the informational columns end) 
Generates an array of kinship matrices based on LOCO, and the corresopnding genotype array by chromosome.
"""


function loco_processing(geno::DataFrame, genome_start_index::Int)
    chr_names = unique(geno[:,"Chr"])
    geno_array = Vector{Matrix{Float64}}(undef, length(chr_names))
    kinship_array = Vector{Matrix{Float64}}(undef, length(chr_names))
    for i in 1:length(chr_names)
        chr = chr_names[i]
        only_chr = Matrix{Float64}(permutedims(subset(geno, :Chr => ByRow(==(chr)))[:,genome_start_index:end]))
        kinship_array[i] = calcKinship(Matrix{Float64}(permutedims(subset(geno, :Chr => ByRow(!=(chr)))[:,genome_start_index:end])))
        geno_array[i] = only_chr
    end
    return geno_array, kinship_array
end


"""
loco_scan(y, geno_array, kinship_array; reml=false, permutation_test=true, nperms=1000, weights=missing, prior_variance=0.0, prior_sample_size=0.0)

Single trait scan without covariates for LOCO data structure.
Takes in the phenotype column matrix y and the geno and kinship arrays that would be returned by the loco_processing function, as well as optional keywords pertaining to the scan.

___

loco_scan(y, geno_array, covar, kinship_array; reml=false, permutation_test=true, nperms=1000, weights=missing, prior_variance=0.0, prior_sample_size=0.0)

Single trait scan with covariates for LOCO data structure.
Takes in the phenotype column matrix y, the covariate column matrix covar, and the geno and kinship arrays that would be returned by the loco_processing function, as well as optional keywords pertaining to the scan.

"""


function loco_scan(y, geno_array, kinship_array; reml=false, permutation_test=true, nperms=1000, weights=missing, prior_variance=0.0, prior_sample_size=0.0)
    results_by_chr = map((g,k) -> scan(y, g, k; reml=reml, permutation_test=permutation_test, nperms=nperms, weights=weights, prior_variance=prior_variance, prior_sample_size=prior_sample_size), geno_array,kinship_array)
    flattened = collect(Iterators.Flatten(results_by_chr))
    results_full = (sigma2_e = reduce(vcat,flattened[1:4:end]), h2_null = reduce(vcat,flattened[2:4:end]), lod = reduce(vcat, flattened[3:4:end]), L_perms = reduce(vcat, flattened[4:4:end]))
    return results_full
end


function loco_scan(y, geno_array, covar, kinship_array; reml=false, permutation_test=true, nperms=1000, weights=missing, prior_variance=0.0, prior_sample_size=0.0)
    results_by_chr = map((g,k) -> scan(y, g, covar, k; reml=reml, permutation_test=permutation_test, nperms=nperms, weights=weights, prior_variance=prior_variance, prior_sample_size=prior_sample_size), geno_array,kinship_array)
    flattened = collect(Iterators.Flatten(results_by_chr))
    results_full = (sigma2_e = reduce(vcat,flattened[1:4:end]), h2_null = reduce(vcat,flattened[2:4:end]), lod = reduce(vcat, flattened[3:4:end]), L_perms = reduce(vcat, flattened[4:4:end]))
    return results_full
end

