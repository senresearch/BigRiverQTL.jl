"""
parse_json(file::String)

Single trait scan without covariates for LOCO data structure.

# Arguments

- `y`is the phenotype column matrix. 
- `dfG`is is a dataframe containing genotype values and genotype info such as the chromosome, loci... 
- `kwargs` are optional keywords arguments pertaining to the `BulkLMM.bulkscan()` function. For
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

"""
function parse_json(file::String)
	indict = open(file, "r") do f
		JSON.parse(f)
	end
	return indict
end


"""

"""
function read_data(filename)
	# read the file into lines
	lines = readlines(filename)
	# which lines have # as first character
	firstpound = (x->match(r"^#",x)).( lines ) .!= nothing
	# last line of comment
	startdata = findfirst(firstpound.==0)

	return CSV.read(filename, DataFrame; header=startdata)
end

