
"""
get_gmap(filename::String)

Creates a `Gmap` type/struct from gmap CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the gmap CSV file.

# Output

Returns a `Gmap` type/struct.


"""
function get_gmap(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	#  check gmap file exists
	if (in("gmap", keys(jsondict)))
		gmapfile = joinpath(data_dir, jsondict["gmap"])
	else
		throw("Error: gmap not found in control file")
	end

	# load gmap file   
	gdf = groupby(read_data(gmapfile), :chr)

	# chromosome
	chr = [group.chr[1] for group in gdf]

	# marker
	marker = [String.(group.marker) for group in gdf]

	# relative position
	pos = [group.pos for group in gdf]

	return Gmap(chr, marker, pos)
end




"""
get_crosstype(file::String)

Creates a `CrossType` type/struct from  control file in json format.

# Argument

- `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `CrossType` type/struct.


"""
function get_crosstype(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   
	jsondict = BigRiverQTL.parse_json(filename)

	# make type

	#crosstype
	if (in("crosstype", keys(jsondict)))
		crosstype = jsondict["crosstype"]
	else
		throw("Error: CrossType not found")
	end

	return CrossType(crosstype)
end


"""
get_alleles(file::String)

Creates a `Alleles` type/struct from  control file in json format.

# Argument

- `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `Alleles` type/struct.
"""
function get_alleles(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   
	jsondict = BigRiverQTL.parse_json(filename)

	# make type

	#alleles
	alleles = jsondict["alleles"]

	return Alleles(alleles)
end


"""
get_genotype(file::String)

Creates a `GenoType` type/struct from  control file in json format.

# Argument

- `filename` : A string containing the name(with directory) of the control file in json format.

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

- `filename` : A string containing the name(with directory) of the control file in json format.

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

- `filename` : A string containing the name(with directory) of the control file in json format.


# Output

Returns a `Geno` type/struct.
"""
function get_geno(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   
	jsondict = BigRiverQTL.parse_json(filename)

	#  check gmap file exists
	if (in("gmap", keys(jsondict)))
		gmapfile = joinpath(data_dir, jsondict["gmap"])
	else
		throw("Error: gmap not found in control file")
	end

	# check geno file exists
	if (in("geno", keys(jsondict)))
		genofile = joinpath(data_dir, jsondict["geno"])
	else
		throw("Error: geno not found in control file")
	end

	# get gmap object
	gmap = get_gmap(filename)

	# load gmap dataframe
	df_gmap = read_data(gmapfile)

	# load geno dataframe
	df_geno = read_data(genofile)

	# get samples names
	samples = names(df_geno)[2:end]

	# get chromosomes names
	chr = gmap.chr

	# get markers names
	marker = gmap.marker_name

	#genotype
	genotype = get_genotype(filename)

	# values
	# get encoded genetotypes
	df_encoded = DataFrame(
		encode_genotype(genotype.label, df_geno[:, 2:end] |> x -> Matrix(x)),
		samples,
	)

	# built encoded geno DataFrames
	insertcols!(df_encoded, 1, :marker => df_geno.marker)

	gdf = innerjoin(
		select(df_gmap, [:chr, :marker]), 
		df_encoded, 
		on = :marker
	) |> x -> groupby(x, :chr);
		
	val = [Matrix(select(group, Not(:marker, :chr))) for group in gdf];
	
	#crosstype
	crosstype = get_crosstype(filename)

	#alleles
	alleles = get_alleles(filename)

	#genotranspose
	genotranspose = get_genotranspose(filename)

	return Geno(samples, chr, marker, val, crosstype, alleles, genotype, genotranspose)
end


"""
get_pmap(filename::String)

Creates a `Pmap` type/struct from Pmap CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the pmap CSV file.

# Output

Returns a `Pmap` type/struct.

"""
function get_pmap(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	#  check pmap file exists
	if (in("pmap", keys(jsondict)))
		pmapfile = joinpath(data_dir, jsondict["pmap"])
	else
		throw("Error: pmap not found in control file")
	end

	# load pmap file   
	gdf = groupby(read_data(pmapfile), :chr)

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

- `filename` : A string containing the name(with directory) of the pheno CSV file.

# Output

Returns a `Pheno` type/struct.
"""
function get_pheno(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load control file   
	jsondict = parse_json(filename)

	#  check pheno file exists
	if (in("pheno", keys(jsondict)))
		phenofile = joinpath(data_dir, jsondict["pheno"])
	else
		throw("Error: pheno not found in control file")
	end	

	# load file   
	df_pheno = read_data(phenofile)

	# samples
	samples = df_pheno[:, 1]

	# traits' names
	traits = names(df_pheno)[2:end]

	# value of the traits
	df_pheno = Matrix(df_pheno)
	val = tryparse.(Float64, df_pheno[:, 2:end])

	return Pheno(samples, traits, val)
end


"""
get_phenocovar(filename::String)

Creates a `Phenocovar` type/struct from phenocovar CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the phenocovar CSV file.

# Output

Returns a `Phenocovar` type/struct.
"""
function get_phenocovar(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   

	df_phenocovar = read_data(filename)

	# make type

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

- `filename` : A string containing the name(with directory) of the crossinfo CSV file.

# Output

Returns a `Crossinfo` type/struct.
"""
function get_crossinfo(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   
	df_crossdirection = read_data(filename)

	# make type

	# traits' names
	samples = df_crossdirection[:, 1]


	# description
	direction = df_crossdirection[:, 2]

	return CrossInfo(samples, direction)
end


"""
get_isxchar(filename::String)

Creates a `IsXChar` type/struct from gmap CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the gmap CSV file.

# Output

Returns a `IsXChar` type/struct.
"""
function get_isxchar(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   

	gmap = get_gmap(filename)

	# make type

	# chromosome
	chr = gmap.chr

	# isXchar

	isxchar = zeros(Bool, length(chr))
	isxchar[findfirst(x -> x == "X", chr)] = 1

	return IsXChar(chr, isxchar)
end


"""
get_isfemale(filename::String)

Creates a `IsFemale` type/struct from control file in json format.

# Argument

- `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `IsFemale` type/struct.
"""
function get_isfemale(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   
	jsondict = BigRiverQTL.parse_json(filename)
	crossinfofile = joinpath(data_dir, jsondict["cross_info"]["file"])

	# make type

	# samples
	samples = get_crossinfo(crossinfofile).sample_id

	# isfemale

	isfemale = ones(Bool, length(samples))
	if (in("sex", keys(jsondict)))
		isfemale[findall(x -> x == "Male", jsondict["sex"])] = 0
	end

	return IsFemale(samples, isfemale)
end


"""
get_geneticstudydata(filename::String)

Creates a `GeneticStudyData` type/struct from  control file in json format.

# Argument

- `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `GeneticStudyData` type/struct.
"""
function get_geneticstudydata(filename::String)
	# check if filename is directory or file name
	data_dir, filename = get_control_file(filename)

	# load file   
	jsondict = BigRiverQTL.parse_json(filename)

	# make type

	#gmap
	if (in("gmap", keys(jsondict)))
		gmapfile = joinpath(data_dir, jsondict["gmap"])

	else
		throw("Error: gmap not found in control file")
	end

	gmap = get_gmap(gmapfile)

	#geno
	if (in("geno", keys(jsondict)))
		genofile = joinpath(data_dir, jsondict["geno"])

	else
		throw("Error: geno not found in control file")
	end

	geno = get_geno(filename)

	#pmap
	if (in("pmap", keys(jsondict)))
		pmapfile = joinpath(data_dir, jsondict["pmap"])

	else
		throw("Error: pmap not found in control file")
	end

	pmap = get_pmap(pmapfile)

	#pheno
	if (in("pheno", keys(jsondict)))
		phenofile = joinpath(data_dir, jsondict["pheno"])

	else
		throw("Error: pheno not found in control file")
	end

	pheno = get_pheno(phenofile)

	#phenocov
	if (in("phenocovar", keys(jsondict)))
		phenocovfile = joinpath(data_dir, jsondict["phenocovar"])

	else
		throw("Error: phenocovar not found in control file")
	end

	phenocov = get_phenocovar(phenocovfile)

	#isXchar
	isXchar = get_isxchar(gmapfile)

	#isfemale
	isfemale = get_isfemale(filename)

	#crossinfo
	crossinfofile = joinpath(data_dir, jsondict["cross_info"]["file"])
	crossinfo = get_crossinfo(crossinfofile)

	return GeneticStudyData(gmap,
		geno,
		pmap,
		pheno,
		phenocov,
		isXchar,
		isfemale,
		crossinfo)
end

