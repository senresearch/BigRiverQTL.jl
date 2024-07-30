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
#     display_name: Julia 1.10.3
#     language: julia
#     name: julia-1.10
# ---

using Pkg
Pkg.activate("../")

using Revise

using BigRiverQTLPlots
using Random, Statistics
using Plots
using Helium
using BigRiverQTL
using CSV
using DataFrames
using BulkLMM

# +
##############
# BXD spleen #
##############

########
# Data #
########
data_dir = joinpath(@__DIR__, "../data/BXD/");
file = joinpath(data_dir, "bxd.json");
# -

# Transforming data to a optimised and accessible data type
data = get_geneticstudydata(file);

gInfo=data.gmap;
pInfo=data.phenocov;
pheno=data.pheno;
pheno=data.pheno.val;
geno=data.geno.val[1];
geno_processed=convert(Array{Float64}, geno);

#################
# Preprocessing #
#################
traitID = 1112;
pheno_y = pheno[:, traitID];
pheno_y2=ones(length(pheno_y));
pheno_y2[findall(x->x!=nothing,pheno_y)]=pheno_y[findall(x->x!=nothing,pheno_y)];

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
########
# Data #
########
bulklmmdir = dirname(pathof(BulkLMM));

gmap_file = joinpath(bulklmmdir, "..", "data", "bxdData", "gmap.csv");
gInfo2 = BulkLMM.CSV.read(gmap_file, BulkLMM.DataFrames.DataFrame);

# -

gInfo2;

# +
function gmap2df(gmap::Gmap)
    Locus=reduce(vcat,gmap.marker_name);
    Mb=reduce(vcat,gmap.pos);
    start=0
    Chr=repeat(["0"],sum(length.(gmap.marker_name)))
    Cm=repeat(["Missing"],sum(length.(gmap.marker_name)))
    for i in 1: eachindex(gmap.chr)
        l_i=length(gmap.marker_name[i])
        Chr[start+1:l_i+start] .= gmap.chr[i]
        start=start+l_i
    end
    df=DataFrame((Locus=Locus,Chr=Chr,Mb=Mb,Cm=Cm))
    return df

end
# -

gmap=gInfo

gmap.chr

Locus=reduce(vcat,gmap.marker_name);
Mb=reduce(vcat,gmap.pos);
start=0
Chr=repeat(["0"],sum(length.(gmap.marker_name)))
Cm=repeat(["Missing"],sum(length.(gmap.marker_name)))
for i in  eachindex(gmap.chr)
    l_i=length(gmap.marker_name[i])
    Chr[start+1:l_i+start] .= gmap.chr[i]
    start=start+l_i
end


gmap2df(gInfo);

# +
function gmap2df(gmap::Gmap)
    Locus=reduce(vcat,gmap.marker_name);
    Mb=reduce(vcat,gmap.pos);
    start=0
    Chr=repeat(["0"],sum(length.(gmap.marker_name)))
    Cm=repeat([Missing],sum(length.(gmap.marker_name)))
    for i in  eachindex(gmap.chr)
        l_i=length(gmap.marker_name[i])
        Chr[start+1:l_i+start] .= gmap.chr[i]
        start=start+l_i
    end
    df=DataFrame(Locus=Locus,Chr=Chr,Mb=Mb,cM=Cm)
    return df

end
# -

gInfo3=gmap2df(gInfo)

# +
########
# Plot #
########

plot_QTL(single_results_perms, gInfo3);


