
"""
get_gmap(filename::String)

Creates a `Gmap` type/struct from gmap CSV file.

# Argument

* `filename` : A string containing the name(with directory) of the gmap CSV file.

# Output

Returns a `Gmap` type/struct.


"""
function get_gmap(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	#  get gmap file
	gmapfile = joinpath(data_dir, check_key(jsondict, "gmap"))

	# load gmap file   
	gdf = groupby(
		read_data(
			gmapfile;
			delim = check_key(jsondict, "sep"),
			comment = check_key(jsondict, "comment.char"),
		),
		:chr,
	)

	# chromosome
	chr = [group.chr[1] for group in gdf]

	# marker
	marker = [String.(group.marker) for group in gdf]

	# relative position
	pos = [group.pos for group in gdf]

		# unit
		unit = ""
		f = open(gmapfile, "r")
	
		s = readline(f)
		s_lower = lowercase(s)  # Convert the input string to lowercase
		if occursin("mbp", s_lower)
			unit = "Mbp"
		elseif occursin("mb", s_lower)
			unit = "Mb"
		elseif occursin("cm", s_lower)
			unit = "cM"
		else
			@warn "No unit detected!"
		end
	
		close(f)

	return Gmap(chr, marker, pos,  unit)
end




"""
get_crosstype(file::String)

Creates a `CrossType` type/struct from  control file in json format.

# Argument

* `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `CrossType` type/struct.


"""
function get_crosstype(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = BigRiverQTL.parse_json(filename)

	# get crosstype
	crosstype = check_key(jsondict, "crosstype")

	return CrossType(crosstype)
end


"""
get_alleles(file::String)

Creates a `Alleles` type/struct from  control file in json format.

# Argument

* `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `Alleles` type/struct.
"""
function get_alleles(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = BigRiverQTL.parse_json(filename)

	#alleles
	alleles = check_key(jsondict, "alleles")

	return Alleles(alleles)
end


"""
get_genotype(file::String)

Creates a `GenoType` type/struct from  control file in json format.

# Argument

* `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `GenoType` type/struct.
"""
function get_genotype(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   
	jsondict = BigRiverQTL.parse_json(filename)

	# genotype
	if (in("genotypes", keys(jsondict)))
		label = jsondict["genotypes"]
	else
		# Rqtl2 default
		label = Dict{String, Int}("A" => 1, "H" => 2, "B" => 3, "D" => 4, "C" => 5)
	end

	return GenoType(label)
end


"""
get_genotranspose(file::String)

Creates a `GenoTranspose` type/struct from  control file in json format.

# Argument

* `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `GenoTranspose` type/struct.
"""
function get_genotranspose(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   
	jsondict = BigRiverQTL.parse_json(filename)

	#g enotype
	if (in("geno_transposed", keys(jsondict)))
		val = jsondict["geno_transposed"]
	else
		val = FALSE
	end

	return GenoTranspose(val)
end


"""
get_geno(filename::String)

Creates a `Geno` type/struct from gmap CSV file and geno CSV file.

# Argument

* `filename` : A string containing the name(with directory) of the control file in json format.


# Output

Returns a `Geno` type/struct.
"""
function get_geno(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	# get gmap file
	gmapfile = joinpath(data_dir, check_key(jsondict, "gmap"))

	# get geno file 
	genofile = joinpath(data_dir, check_key(jsondict, "geno"))

	# get gmap object
	gmap = get_gmap(filename)

	# load gmap dataframe
	df_gmap = read_data(
		gmapfile;
		delim = check_key(jsondict, "sep"),
		comment = check_key(jsondict, "comment.char"),
	)

	# load geno dataframe
	df_geno = read_data(
		genofile;
		delim = check_key(jsondict, "sep"),
		comment = check_key(jsondict, "comment.char"),
	)

	# check if geno data are transposed
	if !jsondict["geno_transposed"]
		# permutedims
		df_geno = permutedims(df_geno, 1, "marker")
	end

	# get samples names
	samples = names(df_geno)[2:end]
	# get values 
	mat_geno = df_geno[:, 2:end] |> x -> Matrix(x) 
	# get markers names
	marker = df_geno[:, 1]

	# # check if geno data are transposed
	# if jsondict["geno_transposed"]
	# 	# get samples names
	# 	samples = names(df_geno)[2:end]
	# 	# get values 
	# 	mat_geno = df_geno[:, 2:end] |> x -> Matrix(x) 
	# 	# get markers names
	# 	marker = df_geno.marker
	# else
	# 	# get samples names
	# 	samples = df_geno[:, 1]
	# 	# get values
	# 	mat_geno = df_geno[:, 2:end] |> x -> Matrix(x) |> permutedims
	# 	# get markers names
	# 	marker = names(df_geno)[2:end]
	# end

	# get chromosomes names
	chr = gmap.chr

	#genotype
	genotype = get_genotype(filename)

	# values
	# get encoded genetotypes
	df_encoded = DataFrame(
		encode_genotype(genotype.label, mat_geno),
		samples,
	)

	# built encoded geno DataFrames
	insertcols!(df_encoded, 1, :marker => marker)

	gdf = innerjoin(
		select(df_gmap, [:chr, :marker]),
		df_encoded,
		on = :marker,
	) |> x -> groupby(x, :chr)

	val = [permutedims(Matrix(select(group, Not(:marker, :chr)))) for group in gdf]

    # marker
	marker = [String.(group.marker) for group in gdf]

	# crosstype
	crosstype = get_crosstype(filename)

	# alleles
	alleles = get_alleles(filename)

	# genotranspose
	genotranspose = get_genotranspose(filename)

	return Geno(samples, chr, marker, val, crosstype, alleles, genotype, genotranspose)
end


"""
get_pmap(filename::String)

Creates a `Pmap` type/struct from Pmap CSV file.

# Argument

* `filename` : A string containing the name(with directory) of the pmap CSV file.

# Output

Returns a `Pmap` type/struct.

"""
function get_pmap(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	# get pmap file
	pmapfile = joinpath(data_dir, check_key(jsondict, "pmap"))

	# load pmap file
	gdf = groupby(
		read_data(
			pmapfile;
			delim = check_key(jsondict, "sep"),
			comment = check_key(jsondict, "comment.char"),
		),
		:chr,
	)

	# chromosome
	chr = [group.chr[1] for group in gdf]

	# marker
	marker = [group.marker for group in gdf]

	# relative position
	pos = [group.pos for group in gdf]

	# unit
	unit = ""
	f = open(pmapfile, "r")

	s = readline(f)
	s_lower = lowercase(s)  # Convert the input string to lowercase
	if occursin("mbp", s_lower)
		unit = "Mbp"
	elseif occursin("mb", s_lower)
		unit = "Mb"
	elseif occursin("cm", s_lower)
		unit = "cM"
	else
		@warn "No unit detected!"
	end

	close(f)

	return Pmap(chr, marker, pos, unit)
end


"""
get_pheno(filename::String)

Creates a `Pheno` type/struct from Pheno CSV file.

# Argument

* `filename` : A string containing the name(with directory) of the pheno CSV file.

# Output

Returns a `Pheno` type/struct.
"""
function get_pheno(filename::String; samples_col_name::Symbol = :id)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	# get pheno file
	phenofile = joinpath(data_dir, check_key(jsondict, "pheno"))

	# load pheno file   
	df_pheno = read_data(
		phenofile;
		delim = check_key(jsondict, "sep"),
		comment = check_key(jsondict, "comment.char"),
		missingstring = "NA",
	)

	# samples
	samples = df_pheno[:, samples_col_name]

	# traits' names
	traits = names(df_pheno)
	setdiff!(traits, [string(samples_col_name)])

	# value of the traits
	val = Matrix(df_pheno[:, traits])

	return Pheno(samples, traits, val)
end


"""
get_phenocovar(filename::String)

Creates a `Phenocovar` type/struct from phenocovar CSV file.

# Argument

* `filename` : A string containing the name(with directory) of the phenocovar CSV file.

# Output

Returns a `Phenocovar` type/struct.
"""
function get_phenocovar(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	# get phenocovar file
	phenocovarfile = joinpath(data_dir, check_key(jsondict, "phenocovar"))

	# load phenocovar file   
	df_phenocovar = read_data(
		phenocovarfile;
		delim = check_key(jsondict, "sep"),
		comment = check_key(jsondict, "comment.char"),
	)

	# traits' names
	traits = string.(df_phenocovar[:, 1])

	# description
	description = df_phenocovar[:, 2]

	return Phenocov(traits, description)
end


"""
get_crossinfo(filename::String)

Creates a `Crossinfo` type/struct from crossinfo CSV file.

# Argument

* `filename` : A string containing the name(with directory) of the crossinfo CSV file.

# Output

Returns a `Crossinfo` type/struct.
"""
function get_crossinfo(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	# get crossdirection file
	crossinfo_dict = check_key(jsondict, "cross_info")
	crossinfofile = joinpath(data_dir, crossinfo_dict["file"])

	# load crossdirection file   
	df_crossdirection = read_data(
		crossinfofile;
		delim = check_key(jsondict, "sep"),
		comment = check_key(jsondict, "comment.char"),
	)

	# samples names
	samples = df_crossdirection[:, 1]

	# description
	direction = df_crossdirection[:, 2]

	return CrossInfo(samples, direction)
end


"""
get_isxchar(filename::String)

Creates a `IsXChar` type/struct from gmap CSV file.

# Argument

* `filename` : A string containing the name(with directory) of the gmap CSV file.

# Output

Returns a `IsXChar` type/struct.
"""
function get_isxchar(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	# load gmap
	gmap = get_gmap(filename)

	# isXchar
	isxchar = map(x -> x == "X", gmap.chr)

	return IsXChar(gmap.chr, isxchar)
end


"""
get_covar(filename::String)

Creates a `Covar` type/struct from gmap CSV file.

# Argument

* `filename` : A string containing the name(with directory) of the gmap CSV file.

# Output

Returns a `Covar` type/struct.
"""
function get_covar(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	# get covar file
	covarfile = check_key(jsondict, "covar")

	if ismissing(covarfile)
		return Covar(covarfile)
	end

	covarfile = joinpath(data_dir)

	# load covar file   
	df_covar = read_data(
		covarfile;
		delim = check_key(jsondict, "sep"),
		comment = check_key(jsondict, "comment.char"),
		missingstring = check_key(jsondict, "na.strings"),
	)

	return Covar(df_covar)
end


"""
get_isfemale(filename::String)

Creates a `IsFemale` type/struct from control file in json format.
If control file does not stipulate sex information, we assume all female 
(*ref: [qtl2](https://github.com/rqtl/qtl2/blob/main/R/read_cross2.R)* line 229).


# Argument

* `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `IsFemale` type/struct.
"""
function get_isfemale(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	# get samples id
	samples = get_crossinfo(filename).sample_id

	# initiate isfemale vector
	isfemale = ones(Bool, length(samples))

	# check for sex information
	if (in("sex", keys(jsondict)))
		sex_dict = jsondict["sex"]
	else
		@warn "No sex information; assuming all female"
		return IsFemale(samples, isfemale)
	end

	if (in("covar", keys(sex_dict)))
		sex_covar = sex_dict["covar"]
		df_covar = get_covar(filename)
		idx_male = findall(x -> x == sex_dict["male"], df_covar[:, sex_covar])
		isfemale[idx_male] = false
	else
		throw("Error: covar-sex not found in control file")
	end

	return IsFemale(samples, isfemale)
end


"""
get_geneticstudydata(filename::String)

Creates a `GeneticStudyData` type/struct from  control file in json format.

# Argument

* `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `GeneticStudyData` type/struct.
"""
function get_geneticstudydata(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   
	jsondict = parse_json(filename)

	# get data info
	gmap = get_gmap(filename)
	geno = get_geno(filename)
	pmap = get_pmap(filename)
	pheno = get_pheno(filename)
	phenocov = get_phenocovar(filename)
	isXchar = get_isxchar(filename)
	isfemale = get_isfemale(filename)
	crossinfo = get_crossinfo(filename)
	covar = get_covar(filename)

	return GeneticStudyData(
		gmap,
		geno,
		pmap,
		pheno,
		phenocov,
		isXchar,
		isfemale,
		crossinfo,
		covar,
	)
end

