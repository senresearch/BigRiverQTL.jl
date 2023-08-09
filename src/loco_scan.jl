"""
loco_scan(y::Matrix{Float64}, dfG::DataFrame; kwargs...)

Single trait scan without covariates for LOCO data structure.

# Arguments

- `y`is the phenotype column matrix. 
- `dfG`is is a dataframe containing genotype values and genotype info such as the chromosome, loci... 
- `kwargs` are optional keywords arguments pertaining to the `BulkLMM.scan()` function. For
	example: `reml` = false, `permutation_test` = true, `nperms` = 1000, 
	`weights` = missing, `prior_variance` = 0.0, `prior_sample_size` = 0.0. Refer to 
	`BulkLMM` documentation for more details.

# Example

```julia
loco_scan(y, df_geno;
		 reml = false, permutation_test = true, nperms = 1000, 
		 weights = missing, prior_variance = 0.0, prior_sample_size = 0.0)
```
___

loco_scan(y::Matrix{Float64}, G::Vector{Matrix{Float64}}, K::Vector{Matrix{Float64}};
	kwargs...)

Single trait scan without covariates for LOCO data structure.

# Arguments

- `y`is the phenotype column matrix. 
- `G`is vector of genotype matrices based on the chromosome.
- `K` is a vector of kinship matrices.
- `kwargs` are optional keywords arguments pertaining to the `BulkLMM.scan()` function. For
	example: `reml` = false, `permutation_test` = true, `nperms` = 1000, 
	`weights` = missing, `prior_variance` = 0.0, `prior_sample_size` = 0.0. Refer to 
	`BulkLMM` documentation for more details.
	

```julia
loco_scan(y, arr_geno, arr_kinship;
		 reml = false, permutation_test = true, nperms = 1000, 
		 weights = missing, prior_variance = 0.0, prior_sample_size = 0.0)
```
___

loco_scan(y::Matrix{Float64}, G::Vector{Matrix{Float64}}, covar::Matrix{Float64},
        K::Vector{Matrix{Float64}};	kwargs...)

Single trait scan with covariates for LOCO data structure.

# Arguments

- `y`is the phenotype column matrix. 
- `G`is vector of genotype matrices based on the chromosome.
- `covar` is covariate column matrix.
- `K` is a vector of kinship matrices.
- `kwargs` are optional keywords arguments pertaining to the `BulkLMM.scan()` function. For
	example: `reml` = false, `permutation_test` = true, `nperms` = 1000, 
	`weights` = missing, `prior_variance` = 0.0, `prior_sample_size` = 0.0. Refer to 
	`BulkLMM` documentation for more details.
	
# Example

```julia
loco_scan(y, arr_geno, covar, arr_kinship;
		 reml = false, permutation_test = true, nperms = 1000, 
		 weights = missing, prior_variance = 0.0, prior_sample_size = 0.0)
```

"""
function loco_scan(y::Matrix{Float64}, df_geno::DataFrame; kwargs...)

	G = get_loco_geno(df_geno; kwargs...)
	N = length(G)
	K = calcLocoKinship(G)

	return loco_scan(y, G, K; kwargs...)
end

function loco_scan(y::Matrix{Float64}, G::Vector{Matrix{Float64}}, K::Vector{Matrix{Float64}};
	kwargs...)

    N = length(G)
	results_loco = map((g, k) -> scan(y, g, k; kwargs...), G, K)

	return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
		h2_null = [results_loco[i].h2_null for i in 1:N],
		lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
		L_perms = reduce(vcat, ([results_loco[i].L_perms for i in 1:N])),
	)
end

function loco_scan(y::Matrix{Float64}, G::Vector{Matrix{Float64}}, covar::Matrix{Float64},
	K::Vector{Matrix{Float64}}; kwargs...)

    N = length(G)
	results_loco = map((g, k) -> scan(y, g, covar, k; kwargs...), G, K)

	return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
		h2_null = [results_loco[i].h2_null for i in 1:N],
		lod = reduce(vcat, ([results_loco[i].lod[2:end] for i in 1:N])),
		L_perms = reduce(vcat, ([results_loco[i].L_perms[2:end, :] for i in 1:N])),
	)
end





function loco_scan2(y, geno_array, kinship_array; reml = false, permutation_test = true, nperms = 1000, weights = missing, prior_variance = 0.0, prior_sample_size = 0.0)
	results_by_chr = map((g, k) -> scan(y, g, k; reml = reml, permutation_test = permutation_test, nperms = nperms, weights = weights, prior_variance = prior_variance, prior_sample_size = prior_sample_size), geno_array, kinship_array)
	flattened = collect(Iterators.Flatten(results_by_chr))
	results_full = (sigma2_e = reduce(vcat, flattened[1:4:end]), h2_null = reduce(vcat, flattened[2:4:end]), lod = reduce(vcat, flattened[3:4:end]), L_perms = reduce(vcat, flattened[4:4:end]))
	return results_full
end


function loco_scan2(y, geno_array, covar, kinship_array; reml = false, permutation_test = true, nperms = 1000, weights = missing, prior_variance = 0.0, prior_sample_size = 0.0)
	results_by_chr = map((g, k) -> scan(y, g, covar, k; reml = reml, permutation_test = permutation_test, nperms = nperms, weights = weights, prior_variance = prior_variance, prior_sample_size = prior_sample_size), geno_array, kinship_array)
	flattened = collect(Iterators.Flatten(results_by_chr))
	results_full = (sigma2_e = reduce(vcat, flattened[1:4:end]), h2_null = reduce(vcat, flattened[2:4:end]), lod = reduce(vcat, flattened[3:4:end]), L_perms = reduce(vcat, flattened[4:4:end]))
	return results_full
end
