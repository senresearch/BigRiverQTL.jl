# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
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


# +
# Data types
# gmap contains 
# makers info 
gInfo = data.gmap;
# pehnotype info 
pInfo = data.phenocov;
# phenotype values 
pheno = data.pheno.val;


# We can get the genotype matrix using the following command:
geno = reduce(hcat, data.geno.val);
# For computing reasons, we need to convert the geno matrix in Float64
geno_processed = geno .- 1.0 |>
    x -> replace(x, missing => 0.5) |>

# geno_processed = convert(Matrix{Float64}, geno_processed);

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
