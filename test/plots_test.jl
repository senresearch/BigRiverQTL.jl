########
data_dir = joinpath(@__DIR__, "data/BXD/");
file = joinpath(data_dir, "bxd.json");


# Transforming data to a optimised and accessible data type
data = get_geneticstudydata(file);

# Data types
gInfo = data.gmap;
pInfo = data.phenocov;
# pheno=data.pheno;
pheno = data.pheno.val;
geno = reduce(vcat, data.geno.val);
geno_processed = geno .- 1.0
replace!(geno_processed, missing => 0.5);
geno_processed = convert(Matrix{Float64}, geno_processed);
geno_processed = permutedims(geno_processed);

#################
# Preprocessing #
#################
traitID = 1112;
pheno_y = pheno[:, traitID];
pheno_y2 = ones(length(pheno_y));
idx_not_missing = findall(!ismissing, pheno_y)
pheno_y2[idx_not_missing] = pheno_y[idx_not_missing];

###########
# Kinship #
###########
kinship = kinship_gs(geno_processed, 0.99);


########
# Scan #
########

single_results_perms = BigRiverQTL.BulkLMM.scan(
	pheno_y2,
	geno_processed,
	kinship;
	permutation_test = true,
	nperms = 1000,
);


#########
# Plots #
#########


# QTL plots
p1 = plot_QTL(single_results_perms, gInfo, mbColname = "Pos");

# Manhattan plots
p2 = plot_manhattan(single_results_perms, gInfo, mbColname = "Pos");

@testset "QTL plot Tests" begin
	@test isa(p1[1][4], Plots.Series)
end

@testset "Mahattan plot Tests" begin
	@test isa(p2[1][2], Plots.Series)
end


