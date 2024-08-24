"""
loco_bulkscan(y::Matrix{Float64}, dfG::DataFrame; kwargs...)

Single trait scan without covariates for LOCO data structure.

# Arguments

* `y`is the phenotype column matrix. 
* `dfG`is is a dataframe containing genotype values and genotype info such as the chromosome, loci... 
* `kwargs` are optional keywords arguments pertaining to the `BulkLMM.bulkscan()` function. For
	example: 
	- method::String = "null-grid", h2_grid::Array{Float64, 1} = collect(0.0:0.1:0.9),
	- nb::Int64 = Threads.nthreads(), 
	- nt_blas::Int64 = 1, 
	- weights::Union{Missing, Array{Float64, 1}} = missing,
	- prior_variance::Float64 = 1.0, prior_sample_size::Float64 = 0.0,
	- reml::Bool = false, optim_interval::Int64 = 1,
	- output_pvals::Bool = false, chisq_df::Int64 = 1 
Refer to `BulkLMM` documentation for more details.

# Example

```julia
loco_bulkscan(y, df_geno;
		 method = "null-grid", h2_grid = collect(0.0:0.1:0.9),
		 reml = false, weights = missing, prior_variance = 0.0, prior_sample_size = 0.0)
```
___

loco_bulkscan(y::Matrix{Float64}, G::Vector{Matrix{Float64}}, K::Vector{Matrix{Float64}};
	kwargs...)

Single trait scan without covariates for LOCO data structure.

# Arguments

* `y`is the phenotype column matrix. 
* `G`is vector of genotype matrices based on the chromosome.
* `K` is a vector of kinship matrices.
* `kwargs` are optional keywords arguments pertaining to the `BulkLMM.scan()` function. For
	example: 
	- method::String = "null-grid", h2_grid::Array{Float64, 1} = collect(0.0:0.1:0.9),
	- nb::Int64 = Threads.nthreads(), 
	- nt_blas::Int64 = 1, 
	- weights::Union{Missing, Array{Float64, 1}} = missing,
	- prior_variance::Float64 = 1.0, prior_sample_size::Float64 = 0.0,
	- reml::Bool = false, optim_interval::Int64 = 1,
	- output_pvals::Bool = false, chisq_df::Int64 = 1 
Refer to `BulkLMM` documentation for more details.
	

```julia
loco_bulkscan(y, arr_geno, arr_kinship;
	method = "null-grid", h2_grid = collect(0.0:0.1:0.9),
	reml = false, weights = missing, prior_variance = 0.0, prior_sample_size = 0.0)
```
___

loco_bulkscan(y::Matrix{Float64}, G::Vector{Matrix{Float64}}, covar::Matrix{Float64},
        K::Vector{Matrix{Float64}};	kwargs...)

Single trait scan with covariates for LOCO data structure.

# Arguments

* `y`is the phenotype column matrix. 
* `G`is vector of genotype matrices based on the chromosome.
* `covar` is covariate column matrix.
* `K` is a vector of kinship matrices.
* `kwargs` are optional keywords arguments pertaining to the `BulkLMM.scan()` function. For
	example: 
	- method::String = "null-grid", h2_grid::Array{Float64, 1} = collect(0.0:0.1:0.9),
	- nb::Int64 = Threads.nthreads(), 
	- nt_blas::Int64 = 1, 
	- weights::Union{Missing, Array{Float64, 1}} = missing,
	- prior_variance::Float64 = 1.0, prior_sample_size::Float64 = 0.0,
	- reml::Bool = false, optim_interval::Int64 = 1,
	- output_pvals::Bool = false, chisq_df::Int64 = 1 
Refer to `BulkLMM` documentation for more details.
	
# Example

```julia
loco_scan(y, arr_geno, covar, arr_kinship;
	method = "null-grid", h2_grid = collect(0.0:0.1:0.9),
	reml = false, weights = missing, prior_variance = 0.0, prior_sample_size = 0.0)
```

"""
function loco_bulkscan(y::Matrix{Float64}, df_geno::DataFrame; kwargs...)

	G = get_loco_geno(df_geno; kwargs...)
	N = length(G)
	K = calcLocoKinship(G)

	return loco_bulkscan(y, G, K; kwargs...)
end

function loco_bulkscan(Y::Matrix{Float64}, G::Vector{Matrix{Float64}}, K::Vector{Matrix{Float64}};
	kwargs...)

    N = length(G)
	results_loco = map((g, k) -> bulkscan(Y, g, k; kwargs...), G, K)

	results_keys = keys(results_loco[1]);

	if :h2_null_list in results_keys
		if :log10Pvals_mat in results_keys
			return (
				L = reduce(vcat, ([results_loco[i].L for i in 1:N])),
				h2_null_list = permutedims(reduce(hcat, ([results_loco[i].h2_null_list for i in 1:N]))),
				log10Pvals_mat = reduce(vcat, ([results_loco[i].log10Pvals_mat for i in 1:N])),
				Chisq_df = results_loco[1].Chisq_df,
			)
		else
			return (
				L = reduce(vcat, ([results_loco[i].L for i in 1:N])),
				h2_null_list = permutedims(reduce(hcat, ([results_loco[i].h2_null_list for i in 1:N]))),
			)
		end
	else
		if :log10Pvals_mat in results_keys
			return (
				L = reduce(vcat, ([results_loco[i].L for i in 1:N])),
				h2_panel = reduce(vcat, ([results_loco[i].h2_panel for i in 1:N])),
				log10Pvals_mat = reduce(vcat, ([results_loco[i].log10Pvals_mat for i in 1:N])),
				Chisq_df = results_loco[1].Chisq_df,
			)
		else
			return (
				L = reduce(vcat, ([results_loco[i].L for i in 1:N])),
				h2_panel = reduce(vcat, ([results_loco[i].h2_panel for i in 1:N])),
			)
		end	
	end

end

function loco_bulkscan(y::Matrix{Float64}, G::Vector{Matrix{Float64}}, covar::Matrix{Float64},
	K::Vector{Matrix{Float64}}; kwargs...)

    N = length(G)
	results_loco = map((g, k) -> scan(y, g, covar, k; kwargs...), G, K)

	results_keys = keys(results_loco[1]);

	if :h2_null_list in results_keys
		if :log10Pvals_mat in results_keys
			return (
				L = reduce(vcat, ([results_loco[i].L for i in 1:N])),
				h2_null_list = permutedims(reduce(hcat, ([results_loco[i].h2_null_list for i in 1:N]))),
				log10Pvals_mat = reduce(vcat, ([results_loco[i].log10Pvals_mat for i in 1:N])),
				Chisq_df = results_loco[1].Chisq_df,
			)
		else
			return (
				L = reduce(vcat, ([results_loco[i].L for i in 1:N])),
				h2_null_list = permutedims(reduce(hcat, ([results_loco[i].h2_null_list for i in 1:N]))),
			)
		end
	else
		if :log10Pvals_mat in results_keys
			return (
				L = reduce(vcat, ([results_loco[i].L for i in 1:N])),
				h2_panel = reduce(vcat, ([results_loco[i].h2_panel for i in 1:N])),
				log10Pvals_mat = reduce(vcat, ([results_loco[i].log10Pvals_mat for i in 1:N])),
				Chisq_df = results_loco[1].Chisq_df,
			)
		else
			return (
				L = reduce(vcat, ([results_loco[i].L for i in 1:N])),
				h2_panel = reduce(vcat, ([results_loco[i].h2_panel for i in 1:N])),
			)
		end	
	end

end
