"""
get_loco_geno(dfG::DataFrame; chromosome_colname::String = "Chr", idx_start::Int = 5)
					   
Returns a vector of genotype matrices based on the chromosome.

# Arguments
* `dfG` is a dataframe containing genotype values and genotype info such as the chromosome, loci... 
* `chromosome_colname` column name containing chromosome information
* `idx_start` indicates the first column index containing the genotype values 
"""
function get_loco_geno(dfG::DataFrame;
	chromosome_colname::String = "Chr",
	idx_start::Int = 5, # TODO is it better to indicate geno column names?
	kwargs...,
)

	gdf = groupby(dfG, chromosome_colname)

	N = length(gdf)

	G = Vector{Matrix{Float64}}(undef, N)

	for i in 1:N
		G[i] = Matrix{Float64}(permutedims(gdf[i][:, idx_start:end]))
	end

	return G
end


"""
get_loco_geno_info(dfG::DataFrame; chromosome_colname::String = "Chr", idx_info = collect(1:4))
					   
Returns a vector of genotype information dataframes based on the chromosome.

# Arguments
* `dfG` is a dataframe containing genotype values and genotype info such as the chromosome, loci... 
* `chromosome_colname` column name containing chromosome information
* `idx_info` indicates columns containing genotype information (e.g., Locus, Mb...) 
"""
function get_loco_geno_info(dfG::DataFrame;
	chromosome_colname = "Chr",
	idx_info = collect(1:4),
	kwargs...,
)

	select(dfG, idx_info) |>
	x -> groupby(x, chromosome_colname) |>
		 x -> DataFrame.(collect(x))
end


"""
calcLocoKinship(G::Vector{Matrix{Float64}})
					   
Calculates kinship matrices leaving out one chromosome out, and returns a vector 
of kinship matrices per chromosome.

# Arguments
* `G` is a vector of genotype matrices
"""
function calcLocoKinship(G::Vector{Matrix{Float64}})

	N = length(G)
	K = Vector{Matrix{Float64}}(undef, N)

	Q = calcKinship.(G)
	KK = sum(Q)

	for i in 1:N
		K[i] = (KK - Q[i]) ./ (N - 1)
	end

	return K
end

