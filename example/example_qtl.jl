# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Julia 1.11.7
#     language: julia
#     name: julia-1.11
# ---

# # Example QTL
# ___

# This notebook gives an example of QTL analysis using `BigRiverQTL.jl`.

# ## Data

# In this example, we will use a dataset available from the `R/qtl2` package. Specifically, we will use the BXD dataset, which is obtained from the [GeneNetwork](https://genenetwork.org/) website.
#
# You can download the BXD genotype data from the following link:
# [Download BXD Genotype Data](https://raw.githubusercontent.com/rqtl/qtl2data/master/BXD/bxd.zip)
#

# ### Example - BXD 

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


# The current version of `BigRiverQTL` does not have imputation functions.
# Remove the  missing data
data = get_data_completecases(data);

# +
# Data types
# gmap contains 
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
# -

# #### Preprocessing

#################
# Preprocessing #
#################
traitID = 1112;
pheno_y = pheno[:, traitID];
pheno_y2 = ones(length(pheno_y));
idx_not_missing = findall(!ismissing, pheno_y)
pheno_y2[idx_not_missing] = pheno_y[idx_not_missing];

# #### Kinship

###########
# Kinship #
###########
kinship = kinship_gs(geno,.99);

# #### Scan

# +
########
# Scan #
########

single_results_perms = scan(
	pheno_y2,
	geno,
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
savefig(joinpath(@__DIR__, "..", "images", "QTL_thrs_example.png"))
# -

# Manhattan plots
plot_manhattan(single_results_perms, gInfo, mbColname = "Pos")
savefig(joinpath(@__DIR__, "..", "images", "manhattan_thrs_example.png"))
