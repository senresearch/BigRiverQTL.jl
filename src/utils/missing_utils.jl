
"""
	get_geno_completecases(geno_missing::Geno) -> Geno

Identify and remove rows or columns from a genotype data structure `Geno` containing missing values, 
prioritizing those with the highest percentage of missing data.

# Arguments
* `geno_missing::Geno`: A genotype data structure containing possibly missing values.

# Returns
* `Geno`: A new genotype data structure with the same fields as the input, but where rows and columns 
with missing data have been systematically removed until no missing data remains.

# Description

The core of the function involves iteratively identifying and removing rows or columns with the highest percentage of missing data. 
This decision is based on whether the maximum percentage of missing data across rows exceeds that across columns. 
Rows or columns are removed until there are no missing values left in the matrix.

After cleaning, the function reconstructs the `Geno` data structure. 
It ensures the data aligns with the expected format by adjusting the dimensions of the data, regrouping by chromosomes, 
and ensuring marker and sample identifiers are correctly aligned.

"""
function get_geno_completecases(geno_missing::Geno)
	geno = deepcopy(geno_missing)

	# get geno values
	mat_geno = reduce(hcat, geno.val)

	# build dataframe
	df_geno = DataFrame(
		mat_geno,
		reduce(vcat, geno.marker_name),
	)
	insertcols!(df_geno, 1, :samples => geno.sample_id)

	# get ismissing matrix
	mat_geno_missing = ismissing.(mat_geno) .* 1

	# get complete cases
	while sum(mat_geno_missing) > 0

		missing_row_wise = sum(mat_geno_missing, dims = 2) .* 100 / (size(mat_geno, 2)) |> vec
		missing_col_wise = sum(mat_geno_missing, dims = 1) .* 100 / (size(mat_geno, 1)) |> vec

		if (maximum(missing_row_wise) > maximum(missing_col_wise))
			idx_missing = getindex(findmax(missing_row_wise), 2)
			deleteat!(df_geno, idx_missing)
		else
			idx_missing = getindex(findmax(missing_col_wise), 2)
			select!(df_geno, Not(idx_missing))
		end

		mat_geno = Matrix(df_geno)
		mat_geno_missing = ismissing.(mat_geno) .* 1
	end

	# permute dataframe to fit Geno type/struct
	df_geno = permutedims(df_geno, 1, "marker")

	# get samples names
	samples = names(df_geno)[2:end]

	# get df_gmap: chromosomes-markers
	chr = deepcopy(geno.marker_name)
	for i in eachindex(geno.chr)
		chr[i][:] .= geno.chr[i]
	end

	df_gmap = DataFrame(
		chr = reduce(vcat, chr),
		marker = reduce(vcat, geno.marker_name),
	)

	# get data per chromosome
	gdf = innerjoin(
		df_gmap,
		df_geno,
		on = :marker,
	) |> x -> groupby(x, :chr)


	# permute back such rows are samples and columns ar markers
	val = [permutedims(Matrix(select(group, Not(:marker, :chr)))) for group in gdf]

	# marker
	marker = [String.(group.marker) for group in gdf]

	# chromosome
	chr = [group.chr[1] for group in gdf]

	return Geno(
		samples,
		chr,
		marker,
		val,
		geno.cross_type,
		geno.alleles,
		geno.geno_type,
		geno.geno_transpose,
	)
end


"""
	summary_missing(geno_missing::Geno; issorted::Bool = false) -> (missing_per_row::DataFrame, missing_per_col::DataFrame)

Calculate the percentage of missing data for each sample and each marker in a genetic dataset.

# Arguments
* `geno_missing::Geno`: A genetic data structure of type Geno.
* `issorted::Bool = false`: An optional boolean flag indicating whether to sort the results by 
   the percentage of missing data in descending order.

# Returns
* `(missing_per_row::DataFrame, missing_per_col::DataFrame)`: A tuple of two `DataFrame`s:
  * `missing_per_row`: Each row represents a sample with columns for sample identifiers and the percentage of missing data.
  * `missing_per_col`: Each row represents a marker with columns for marker names and the percentage of missing data.
"""
function summary_missing(geno_missing::Geno; issorted::Bool = false)
	# get geno values
	mat_geno = reduce(hcat, geno_missing.val)

	# get ismissing matrix
	mat_geno_missing = ismissing.(mat_geno) .* 1

	missing_row_wise = sum(mat_geno_missing, dims = 2) .* 100 / (size(mat_geno, 2)) |> vec
	missing_col_wise = sum(mat_geno_missing, dims = 1) .* 100 / (size(mat_geno, 1)) |> vec

	# build dataframe
	df_row_missing = DataFrame(
		samples = geno_missing.sample_id,
		percentage = round.(missing_row_wise, digits = 2),
	)

	df_col_missing = DataFrame(
		markers = reduce(vcat, geno_missing.marker_name),
		percentage = round.(missing_col_wise, digits = 2),
	)

	if issorted
		sort!(df_row_missing, [:percentage], rev = true)
		sort!(df_col_missing, [:percentage], rev = true)
	end

	return (missing_per_row = df_row_missing, missing_per_col = df_col_missing)
end