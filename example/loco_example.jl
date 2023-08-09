using BulkLMM
using DataFrames
using BenchmarkTools
using BigRiverQTL

########
# Data #
########
bulklmmdir = dirname(pathof(BulkLMM));

gmap_file = joinpath(bulklmmdir, "..", "data", "bxdData", "gmap.csv");
gInfo = BulkLMM.CSV.read(gmap_file, BulkLMM.DataFrames.DataFrame);

phenocovar_file = joinpath(bulklmmdir, "..", "data", "bxdData", "phenocovar.csv");
pInfo = BulkLMM.CSV.read(phenocovar_file, BulkLMM.DataFrames.DataFrame);

pheno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-pheno-nomissing.csv");
pheno = BulkLMM.DelimitedFiles.readdlm(pheno_file, ',', header = false);
pheno_processed = pheno[2:end, 2:(end-1)] .* 1.0; # exclude the header, the first (transcript ID)and the last columns (sex)

geno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-bxd-genoprob.csv")
geno = BulkLMM.DelimitedFiles.readdlm(geno_file, ',', header = false);
geno_processed = geno[2:end, 1:2:end] .* 1.0;


dfgeno = copy(gInfo)

df_tmp = DataFrame(permutedims(geno_processed), :auto);
insertcols!(df_tmp, 1, :Locus => gInfo.Locus)
df_geno = leftjoin(gInfo, df_tmp, on = :Locus);

############
# Function #
############

vG = get_loco_geno(df_geno);
N = length(vG)
vK = calcLocoKinship(vG);
traitID = 1112;
y = reshape(pheno_processed[:, traitID], :, 1);
@btime begin
results_by_chr = map((g, k) -> scan(y, g, k;
		reml = false,
		permutation_test = true,
		nperms = 1000, weights = missing,
		prior_variance = 0.0,
		prior_sample_size = 0.0),
vG, vK);
sigma2_e = [results_by_chr[i].sigma2_e for i in 1:N]
h2_null = [results_by_chr[i].h2_null for i in 1:N]
lod = reduce(vcat,([results_by_chr[i].lod for i in 1:N]))
L_perms = reduce(vcat,([results_by_chr[i].L_perms for i in 1:N]));
end;

@btime begin
    results_by_chr = map((g, k) -> scan(y, g, k;
            reml = false,
            permutation_test = true,
            nperms = 1000, weights = missing,
            prior_variance = 0.0,
            prior_sample_size = 0.0),
    vG, vK);
    flattened = collect(Iterators.Flatten(results_by_chr))
    sigma2_e = reduce(vcat,flattened[1:4:end])
    h2_null = reduce(vcat,flattened[2:4:end])
    lod = reduce(vcat, flattened[3:4:end])
    L_perms = reduce(vcat, flattened[4:4:end]);
end;

results_full = (sigma2_e = reduce(vcat,flattened[1:4:end]), h2_null = reduce(vcat,flattened[2:4:end]), lod = reduce(vcat, flattened[3:4:end]), L_perms = reduce(vcat, flattened[4:4:end]))


@btime begin
	gd = groupby(df_geno, :Chr)
	N = length(gd)

	G = Vector{Matrix{Float64}}(undef, N)
	K = Vector{Matrix{Float64}}(undef, N)
	for i in 1:N
		G[i] = Matrix{Float64}(permutedims(gd[i][:, 5:end]))
	end

	Q = calcKinship.(geno_array)
	KK = sum(Q)
	for i in 1:N
		K[i] = (KK - Q[i]) ./ (N - 1)
	end
end

genome_start_index = 5
select(df_geno, collect(1:4)) |>
x -> groupby(x, :Chr) |>
	 x -> DataFrame.(collect(x))

@btime begin
	chr_names = unique(df_geno[:, "Chr"])
	geno_array = Vector{Matrix{Float64}}(undef, length(chr_names))
	kinship_array = Vector{Matrix{Float64}}(undef, length(chr_names))
	for i in 1:length(chr_names)
		chr = chr_names[i]
		only_chr = Matrix{Float64}(permutedims(subset(df_geno, :Chr => ByRow(==(chr)))[:, genome_start_index:end]))
		kinship_array[i] = calcKinship(Matrix{Float64}(permutedims(subset(df_geno, :Chr => ByRow(!=(chr)))[:, genome_start_index:end])))
		geno_array[i] = only_chr
	end
end

