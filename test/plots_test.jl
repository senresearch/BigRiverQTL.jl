########
data_dir = joinpath(@__DIR__, "../data/BXD/");
file = joinpath(data_dir, "bxd.json");
# -

# Transforming data to a optimised and accessible data type
data = get_geneticstudydata(file);

# Data types
gInfo=data.gmap;
pInfo=data.phenocov;
pheno=data.pheno;
pheno=data.pheno.val;
geno=reduce(hcat, data.geno.val);
geno_processed=convert(Array{Float64}, geno);

#################
# Preprocessing #
#################
traitID = 1112;
pheno_y = pheno[:, traitID];
pheno_y2=ones(length(pheno_y));
idx_nothing = findall(x->x!=nothing,pheno_y)
pheno_y2[idx_nothing]=pheno_y[idx_nothing];

###########
# Kinship #
###########
kinship = kinship_gs(geno_processed,.99);

# +
########
# Scan #
########

single_results_perms = scan(
	pheno_y2,
	geno_processed,
	kinship;
	permutation_test = true,
	nperms = 1000,
);

# +
#########
# Plots #
#########


# QTL plots
p1=plot_QTL(single_results_perms, gInfo, mbColname = "Pos")




# -

# Manhattan plots
p2=plot_manhattan(single_results_perms, gInfo, mbColname = "Pos")



@testset "QTL plot Tests" begin
    @test isa(p1[1][4],Series)
end

@testset "Mahattan plot Tests" begin
    @test isa(p2[1][3],Series)
end


