
###########################
# Test Plotting functions #
###########################
@testset "Testing plotting function" begin
	########
	data_dir = joinpath(@__DIR__, "data/BXD/");
	file = joinpath(data_dir, "bxd.json");

	# Transforming data to a optimised and accessible data type
	data = get_geneticstudydata(file);

	# Remove the  missing data
	data = get_data_completecases(data);

	# Data types
	# makers info 
	gInfo = data.gmap;

	# phenotype info 
	pInfo = data.phenocov;
	# phenotype values 
	pheno = data.pheno.val;

	# We can get the genotype matrix using the following command.
	# For computing reasons, we need to convert the geno matrix in Float64.
	# One way to do it is to multiply by 1.0
	geno = reduce(hcat, data.geno.val).*1.0;

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
	kinship = kinship_gs(geno, 0.99);


	########
	# Scan #
	########

	single_results_perms = BigRiverQTL.scan(
		pheno_y2,
		geno,
		kinship;
		permutation_test = true,
		nperms = 1000,
	);

	#########
	# Plots #
	#########

	# QTL plots
	p1a = plot_QTL(single_results_perms, gInfo, mbColname = "Pos");
	p1b = plot_QTL(single_results_perms.lod, gInfo, mbColname = "Pos");

	# Manhattan plots
	p2a = plot_manhattan(single_results_perms, gInfo, mbColname = "Pos");
	p2b = plot_manhattan(single_results_perms.lod, gInfo, mbColname = "Pos");

	@testset "QTL plot Tests" begin
		@test isa(p1a[1][4], Plots.Series)
		@test isa(p1b[1][3], Plots.Series)
	end

	@testset "Manhattan plot Tests" begin
		@test isa(p2a[1][2], Plots.Series)
		@test isa(p2b[1][1], Plots.Series)
	end

end
