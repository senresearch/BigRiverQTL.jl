
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

# Transforming data to a optimised and accessible data type
data = get_geneticstudydata(file);

df_g = gmap2df(data.gmap);

df_p1 = pmap2df(data.pmap);

df_p2 = pmap2df(data.pmap);



jsondict = BigRiverQTL.parse_json(file)

# function get_geno(filename::String)
	filename = file
    # check if filename is directory or file name
    data_dir, filename = get_control_file(filename)

    # load control file   
    jsondict = BigRiverQTL.parse_json(filename)

    # get gmap file
    gmapfile = joinpath(data_dir, BigRiverQTL.check_key(jsondict, "gmap"))

    # get geno file 
    genofile = joinpath(data_dir, BigRiverQTL.check_key(jsondict, "geno"))

    # get gmap object
    gmap = get_gmap(filename)

    # load gmap dataframe
    df_gmap = BigRiverQTL.read_data(
        gmapfile;
        delim= BigRiverQTL.check_key(jsondict, "sep"),
        comment=BigRiverQTL.check_key(jsondict, "comment.char")
    )

    # load geno dataframe
    df_geno = BigRiverQTL.read_data(
        genofile;
        delim=BigRiverQTL.check_key(jsondict, "sep"),
        comment=BigRiverQTL.check_key(jsondict, "comment.char")
    )

    # check if geno data are transposed
    if jsondict["geno_transposed"]
        # get samples names
        samples = names(df_geno)[2:end]
        # get values 
        mat_geno = df_geno[:, 2:end] |> x -> Matrix(x)
        # get markers names
        marker = df_geno.marker[:]
    else
        # get samples names
        samples = df_geno[:, 1]
        # get values
        mat_geno = df_geno[:, 2:end] |> x -> Matrix(x) |> permutedims
        # get markers names
        marker =  names(df_geno)[2:end]
    end

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
        on=:marker,
    ) |> x -> groupby(x, :chr)

    val = [Matrix(select(group, Not(:marker, :chr))) for group in gdf]
	# marker
	marker = [String.(group.marker) for group in gdf]
    # crosstype
    crosstype = get_crosstype(filename)

    # alleles
    alleles = get_alleles(filename)

    # genotranspose
    genotranspose = get_genotranspose(filename)

    Geno(samples, chr, marker, val, crosstype, alleles, genotype, genotranspose)


	struct Geno2
		sample_id::Vector{AbstractString}
		chr::Vector{AbstractString}
		marker_name::Vector{Vector{AbstractString}}
	end


	struct Geno4
		sample_id::Vector{AbstractString}
		chr::Vector{AbstractString}
		# marker_name::Vector{Vector{AbstractString}}
	end



	Geno2(samples, chr, marker)
	Geno4(samples, chr)

	gmap.marker