
# https://github.com/rqtl/qtl2data/blob/main/DO_Recla/recla.json

using BigRiverQTL, LinearAlgebra
using Statistics
using Test
using DelimitedFiles
using CSV
using DataFrames
using Plots



data_dir = joinpath(@__DIR__, "data/BXD/");
file = joinpath(data_dir, "bxd.json");
jsondict = BigRiverQTL.parse_json(file)


	#  check pheno file exists
	if (in("pheno", keys(jsondict)))
		phenofile = joinpath(data_dir, jsondict["pheno"])
	else
		throw("Error: pheno not found in control file")
	end	

	# load file 
	
	df_pheno = CSV.read(phenofile, DataFrame; comment="#")

	df_pheno = 	BigRiverQTL.read_data(
		phenofile; 
		missingstring="NA");

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

	unit


# function get_geno(filename::String)
filename = file
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
	gmap = get_gmap(file)

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

	# return Geno(samples, chr, marker, val, crosstype, alleles, genotype, genotranspose)
# end






































# Transforming data to a optimised and accessible data type
# data = get_geneticstudydata(file);

# check if filename is directory or file name
data_dir, filename = BigRiverQTL.get_control_file(file)

# load file   
jsondict = BigRiverQTL.parse_json(filename)


crossinfofile = joinpath(data_dir, jsondict["cross_info"]["file"])

	# make type

	# samples
	df_crossdirection = BigRiverQTL.read_data(crossinfofile)
	get_crossinfo(crossinfofile)


	samples = get_crossinfo(crossinfofile).sample_id

	# isfemale

	isfemale = ones(Bool, length(samples))
	if (in("sex", keys(jsondict)))
		isfemale[findall(x -> x == "Male", jsondict["sex"])] = 0
	end







