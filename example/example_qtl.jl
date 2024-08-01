# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# # Example QTL
# ___

# This notebook gives an example......

# ## Data

# We now want to get genotype data for the BXD panel. We first need to install the R/qtl2 package. As with R/GNapi, it is not available on CRAN, but rather is distributed via a private repository.
#
# [https://raw.githubusercontent.com/rqtl/qtl2data/master/BXD/bxd.zip](https://raw.githubusercontent.com/rqtl/qtl2data/master/BXD/bxd.zip)

# ### Example BXD 

# Libraries
using BigRiverQTL
using Plots

# #### Data

# We assume that the data is stored in `..\data\BXD` directory.

########
# Data #
########
data_dir = joinpath(@__DIR__, "../data/BXD/");
file = joinpath(data_dir, "bxd.json");

# Load bxd data using the function `get_geneticstudydata()`: 

# Load bxd data
data = get_geneticstudydata(file);


# +
# Data types
# gmap contains 
gInfo = data.gmap;
# gmap contains 
pInfo = data.phenocov;
# gmap contains 
pheno = data.pheno;
# gmap contains 
pheno = data.pheno.val;

# We can get the genotype matrix using the following command:
geno = reduce(hcat, data.geno.val);
# For computing reasons, we need to convert the geno matrix in Float64
geno_processed = convert(Array{Float64}, geno);
# -

# #### Preprocessing

#################
# Preprocessing #
#################
traitID = 1112;
pheno_y = pheno[:, traitID];
pheno_y2=ones(length(pheno_y));
idx_nothing = findall(x->x!=nothing,pheno_y)
pheno_y2[idx_nothing]=pheno_y[idx_nothing];

# #### Kinship

###########
# Kinship #
###########
kinship = kinship_gs(geno_processed,.99);

# #### Scan

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
# -

# #### Preprocessing

# #### Plots

# +
#########
# Plots #
#########

# QTL plots
plot_QTL(single_results_perms, gInfo, mbColname = "Pos")

# -

# Manhattan plots
plot_manhattan(single_results_perms, gInfo, mbColname = "Pos")
