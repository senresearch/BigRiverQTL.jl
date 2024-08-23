"""
    subset_geno(
        geno_selection::Geno,
        subset_selection::Union{Vector{String}, InvertedIndex{Vector{String}}};
        isselection_marker::Bool = true
    ) -> Geno

Subset genetic data based on specified markers or samples.

# Arguments
* `geno_selection`: The original genetic data structured as a `Geno` type, which should include properties like sample IDs, 
    marker names, chromosome numbers, etc.
* `subset_selection`: A list of marker or sample 
    names to select, or an `InvertedIndex` for deselection.
* `isselection_marker`: A flag to indicate whether the selection is based on 
    markers (`true`) or samples (`false`). Defaults to `true`.

# Returns
- A new `Geno` object containing the subset of data based on the provided criteria.
"""
function subset_geno(
	geno_selection::Geno,
	subset_selection::Union{Vector{String}, InvertedIndex{Vector{String}}};
	isselection_marker::Bool = true,
)

	geno = deepcopy(geno_selection)

	# get geno values
	mat_geno = reduce(hcat, geno.val)

	# build dataframe
	df_geno = DataFrame(
		mat_geno,
		reduce(vcat, geno.marker_name),
	)
	insertcols!(df_geno, 1, :sample => geno.sample_id)

    # selection
	if isselection_marker
		# selection
		if (typeof(subset_selection) <: InvertedIndex{Vector{String}})
			select!(df_geno, subset_selection)
		else
			select!(df_geno, vcat(["sample"], subset_selection))
		end

		# permute dataframe to fit Geno type/struct
		df_geno = permutedims(df_geno, 1, "marker")
	else
		# permute dataframe to fit Geno type/struct
		df_geno = permutedims(df_geno, 1, "marker")

		# selection
		if (typeof(subset_selection) <: InvertedIndex{Vector{String}})
			select!(df_geno, subset_selection)
		else
			select!(df_geno, vcat(["marker"], subset_selection))
		end
	end

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
	markers = [String.(group.marker) for group in gdf]

    # chromosome
	chr = [group.chr[1] for group in gdf]

	return Geno(
		samples,
		chr,
		markers,
		val,
		geno.cross_type,
		geno.alleles,
		geno.geno_type,
		geno.geno_transpose,
	)
end


"""
    subset_gmap(
        gmap_selection::Gmap,
        subset_selection::Union{Vector{String}, InvertedIndex{Vector{String}}}
    ) -> Gmap

Create a subset of a genetic map (`Gmap`) based on specified markers.

# Arguments
* `gmap_selection`: A genetic map from which a subset will be created. 
It is assumed that `Gmap` is a custom type or struct representing genetic maps, 
which might contain fields such as chromosome (Chr), marker (Locus), and position (Pos).
* `subset_selection`: A list of marker names to select, or an 
`InvertedIndex` for deselection.

# Returns
* A new `Gmap` object representing the subset of the original genetic map. 
The subset contains only the specified markers or excludes specific markers, 
depending on the `subset_selection` parameter.

"""
function subset_gmap(
	gmap_selection::Gmap,
	subset_selection::Union{Vector{String}, InvertedIndex{Vector{String}}}
)
    gmap = deepcopy(gmap_selection)
  
	# build pheno dataframe
	df_gmap = gmap2df(gmap)
    
	# permute dataframe to fit Gmap type/struct
	df_gmap = permutedims(df_gmap, 1, "variables")

    # selection
    if (typeof(subset_selection) <: InvertedIndex{Vector{String}})
        select!(df_gmap, subset_selection)
    else
        select!(df_gmap, vcat(["variables"], subset_selection))
    end

    # permute dataframe to fit Gmap type/struct
	df_gmap = permutedims(df_gmap, 1, "Locus")

    # load gmap file   
	gdf = groupby(df_gmap, :Chr)

	# chromosome
	chr = [group.Chr[1] for group in gdf]

	# marker
	marker = [String.(group.Locus) for group in gdf]

	# relative position
	pos = [group.Pos for group in gdf]

	return Gmap(chr, marker, pos)
end


"""
    subset_pheno(
        pheno_selection::Pheno,
        subset_selection::Union{Vector{String}, InvertedIndex{Vector{String}}}
    ) -> Pheno

Filter and return a subset of a phenotype data structure based on specified criteria.

# Arguments
* `pheno_selection`: A `Pheno` object containing phenotype data.
* `subset_selection`: This argument can either be a vector of strings specifying 
the sample id to retain, or an `InvertedIndex` of a vector of strings 
specifying the sample id to exclude.

# Returns
* A new `Pheno` object containing only the selected subset of traits and 
corresponding samples.

"""
function subset_pheno(
	pheno_selection::Pheno,
	subset_selection::Union{Vector{String}, InvertedIndex{Vector{String}}}
)
    pheno = deepcopy(pheno_selection)
  
	# build pheno dataframe
	df_pheno = DataFrame(
		permutedims(pheno.val),
        pheno.sample_id
	)
	insertcols!(df_pheno, 1, :trait => pheno.traits)
    
    # selection
    if (typeof(subset_selection) <: InvertedIndex{Vector{String}})
        select!(df_pheno, subset_selection)
    else
        select!(df_pheno, vcat(["trait"], subset_selection))
    end

    return Pheno(
		names(df_pheno)[2:end],
		df_pheno.trait[:],
		permutedims(Matrix(df_pheno[:, 2:end]))
	)
end


"""
    subset_isXchar(isXchar_selection::IsXChar, subset_selection::Union{Vector{String}, InvertedIndex{Vector{String}}}) -> IsFemale

Create a subset of an `IsXChar` selection based on specified columns and return an `IsFemale` object containing the subsetted data.

# Arguments
* `isXchar_selection`: An `IsXChar` object containing a data representation which 
presumably includes values (`val`) and chromosome identifiers (`chr`).
* `subset_selection`: This argument can either be a vector of strings specifying 
the chromosomes to retain, or an `InvertedIndex` of a vector of strings 
specifying the chromosomes to exclude.

# Returns
* A new `IsXChar` object constructed with the names of the subsetted dataframe 
as chromosome identifiers and the first row of the subsetted dataframe as values.

"""
function subset_isXchar(
	isXchar_selection::IsXChar,
	subset_selection::Union{Vector{String}, InvertedIndex{Vector{String}}}
)
    isXchar = deepcopy(isXchar_selection)
  
  
	# build isXchar dataframe
	df_isxchar = DataFrame(
        permutedims(isXchar.val),
        isXchar.chr
	)
	   
    # selection
    select!(df_isxchar, subset_selection)
    
    return IsXChar(
		names(df_isxchar),
		Vector(df_isxchar[1, :])
	)
end


function subset_isfemale(
	isfemale_selection::IsFemale,
	subset_selection::Union{Vector{String}, InvertedIndex{Vector{String}}}
)
    isfemale = deepcopy(isfemale_selection)
  
	# build isfemale dataframe
	df_isfemale = DataFrame(
		permutedims(isfemale.val),
        isfemale.sample_id
	)
	   
    # selection
    select!(df_isfemale, subset_selection)
    
    return IsFemale(
		names(df_isfemale),
		Vector(df_isfemale[1, :])
	)
end


function subset_covar(
	covar_selection::Covar,
	subset_selection::Union{Vector{String}, InvertedIndex{Vector{String}}}
)
    # check if covar is missing
    if ismissing(covar_selection.val)
        return Covar(covar_selection.val)
    end

    covar = deepcopy(covar_selection)
  
	# covar dataframe
	df_covar = covar.val

    # permute dataframe to fit Geno type/struct
	df_covar = permutedims(df_covar, 1, "covariates")

    # selection
    if (typeof(subset_selection) <: InvertedIndex{Vector{String}})
        select!(df_covar, subset_selection)
    else
        select!(df_covar, vcat(["covariates"], subset_selection))
    end

    return Covar(permutedims(df_covar, 1, "id"))
end


"""
	select_sample(geno_selection::Geno, markers_selection:: Union{Vector{String}, InvertedIndex{Vector{String}}}) -> Geno

Select and return a subset of genetic samples from a `Geno` object based on specified samples names or their indices.

# Arguments
* `geno_selection`: A `Geno` type or struct that contains genetic data.
* `samples_selection`: A vector of samples names or an inverted index of samples names 
	to be selected from the `geno_selection`. The selection is applied to the samples names in the `Geno` object.

# Returns
* `Geno`: A new `Geno` object that contains only the selected samples, their corresponding genetic information, 
and the subset of the original genotypes related to these samples.
___



"""
function select_sample(geno_selection::Geno, sample_selection::Union{Vector{String}, InvertedIndex{Vector{String}}})
    return subset_geno(geno_selection, sample_selection; isselection_marker = false)
end

function select_sample(gsd_selection::GeneticStudyData, sample_selection::Union{Vector{String}, InvertedIndex{Vector{String}}})
        
    ########
    # Geno #
    ########
    # select sample in geno
    geno_sub = subset_geno(gsd_selection.geno, sample_selection; isselection_marker = false)

    #########
    # Pheno #
    #########
    # select sample in pheno
    pheno_sub = subset_pheno(gsd_selection.pheno, sample_selection)

    ############
    # IsFemale #
    ############
    # select sample in isfemale
    isfemale_sub = subset_isfemale(gsd_selection.isfemale, sample_selection)

    #########
    # Covar #
    #########
    # select sample in covar
    if !ismissing(gsd_selection.covar)
        covar_sub = subset_covar(gsd_selection.covar, sample_selection)
    else
        covar_sub = missing
    end

    return GeneticStudyData(
        gsd_selection.gmap,
        geno_sub,
        gsd_selection.pmap,
        pheno_sub,
        gsd_selection.phenocov,
        gsd_selection.isXchar,
        isfemale_sub,
        gsd_selection.crossinfo,
        covar_sub
    )
end


"""
	select_marker(geno_selection::Geno, markers_selection:: Union{Vector{String}, InvertedIndex{Vector{String}}}) -> Geno

Select and return a subset of genetic markers from a `Geno` object based on specified marker names or their indices.

# Arguments
* `geno_selection`: A `Geno` type or struct that contains genetic data.
* `markers_selection`: A vector of marker names or an inverted index of marker names 
	to be selected from the `geno_selection`. The selection is applied to the marker names in the `Geno` object.

# Returns
* `Geno`: A new `Geno` object that contains only the selected markers, their corresponding genetic information, 
and the subset of the original genotypes related to these markers.

"""
function select_marker(geno_selection::Geno, marker_selection::Union{Vector{String}, InvertedIndex{Vector{String}}})
    return subset_geno(geno_selection, marker_selection; isselection_marker = true)
end

function select_marker(gsd_selection::GeneticStudyData, marker_selection::Union{Vector{String}, InvertedIndex{Vector{String}}})
            
    ########
    # Gmap #
    ########
    # select marker in gmap
    gmap_sub = subset_gmap(gsd_selection.gmap, marker_selection);
    
    ########
    # Geno #
    ########
    # select sample in geno
    geno_sub = subset_geno(gsd_selection.geno, marker_selection; isselection_marker = true);

    ############
    # IsXChar #
    ############
    # select sample in isfemale
    isXchar_sub = subset_isfemale(gsd_selection.isXchar, gmap_sub.chr)

    
    return GeneticStudyData(
        gmap_sub,
        geno_sub,
        gsd_selection.pmap,
        gsd_selection.pheno,
        gsd_selection.phenocov,
        isXchar_sub,
        gsd_selection.isfemale,
        gsd_selection.crossinfo,
        gsd_selection.covar
    )
end


