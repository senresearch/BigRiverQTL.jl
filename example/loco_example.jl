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

G = get_loco_geno(df_geno);
K = calcLocoKinship(G);



@btime begin
gd = groupby(df_geno, :Chr); 
N = length(gd)

G = Vector{Matrix{Float64}}(undef, N)
K = Vector{Matrix{Float64}}(undef, N)
for i in 1:N
    G[i] = Matrix{Float64}(permutedims(gd[i][:,5:end]))
end

Q = calcKinship.(geno_array);
KK = sum(Q)
for i in 1:N
    K[i] = (KK-Q[i])./(N-1)
end
end

genome_start_index = 5
select(df_geno, collect(1:4)) |> 
    x -> groupby(x, :Chr) |>
    x -> DataFrame.(collect(x))

@btime begin
chr_names = unique(df_geno[:,"Chr"])
geno_array = Vector{Matrix{Float64}}(undef, length(chr_names))
kinship_array = Vector{Matrix{Float64}}(undef, length(chr_names))
for i in 1:length(chr_names)
    chr = chr_names[i]
    only_chr = Matrix{Float64}(permutedims(subset(df_geno, :Chr => ByRow(==(chr)))[:,genome_start_index:end]))
    kinship_array[i] = calcKinship(Matrix{Float64}(permutedims(subset(df_geno, :Chr => ByRow(!=(chr)))[:,genome_start_index:end])))
    geno_array[i] = only_chr
end
end

