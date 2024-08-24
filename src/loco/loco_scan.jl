"""
loco_scan(y::Matrix{Float64}, dfG::DataFrame; kwargs...)

Single trait scan without covariates for LOCO data structure.

# Arguments

* `y`is the phenotype column matrix. 
* `dfG`is is a dataframe containing genotype values and genotype info such as the chromosome, loci... 
* `kwargs` are optional keywords arguments pertaining to the `BulkLMM.scan()` function. For
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

* `y`is the phenotype column matrix. 
* `G`is vector of genotype matrices based on the chromosome.
* `K` is a vector of kinship matrices.
* `kwargs` are optional keywords arguments pertaining to the `BulkLMM.scan()` function. For
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

* `y`is the phenotype column matrix. 
* `G`is vector of genotype matrices based on the chromosome.
* `covar` is covariate column matrix.
* `K` is a vector of kinship matrices.
* `kwargs` are optional keywords arguments pertaining to the `BulkLMM.scan()` function. For
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

	results_keys = keys(results_loco[1]);

	if :L_perms in results_keys
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
			L_perms = reduce(vcat, ([results_loco[i].L_perms for i in 1:N])),
		)
	elseif :h2_each_marker in results_keys
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			h2_each_marker = reduce(vcat, ([results_loco[i].h2_each_marker for i in 1:N])),
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
		)
	else
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
		)		
	end
end

function loco_scan(y::Matrix{Float64}, G::Vector{Matrix{Float64}}, covar::Matrix{Float64},
	K::Vector{Matrix{Float64}}; kwargs...)

    N = length(G)
	results_loco = map((g, k) -> scan(y, g, covar, k; kwargs...), G, K)

	results_keys = keys(results_loco[1]);

	if :L_perms in results_keys
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
			L_perms = reduce(vcat, ([results_loco[i].L_perms for i in 1:N])),
		)
	elseif :h2_each_marker in results_keys
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			h2_each_marker = reduce(vcat, ([results_loco[i].h2_each_marker for i in 1:N])),
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
		)
	else
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
		)		
	end
end
