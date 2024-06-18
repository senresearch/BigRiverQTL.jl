# -*- coding: utf-8 -*-
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
#     display_name: Julia_6_Threads 1.10.3
#     language: julia
#     name: julia_6_threads-1.10
# ---

# +

using Statistics
using Distributed
using LinearAlgebra
import StatsBase: sample

# include("Miscellanea.jl")

# -

#genoprob
gpr=[4.44965e-9  1.21848e-6  1.21848e-6  0.999998    8.34238e-16  2.23014e-12  2.23014e-12  1.0          8.34192e-16  2.25839e-12  2.25839e-12  1.0
 1.21848e-6  4.44965e-9  0.999998    1.21848e-6  2.23014e-12  8.34238e-16  1.0          2.23014e-12  2.25839e-12  8.34192e-16  1.0          2.25839e-12
 4.44965e-9  1.21848e-6  1.21848e-6  0.999998    8.34238e-16  2.23014e-12  2.23014e-12  1.0          8.34192e-16  2.25839e-12  2.25839e-12  1.0
 1.21848e-6  0.999998    4.44965e-9  1.21848e-6  2.23014e-12  1.0          8.34238e-16  2.23014e-12  2.25839e-12  1.0          8.34192e-16  2.25839e-12
 0.999998    1.21848e-6  1.21848e-6  4.44965e-9  1.0          2.23014e-12  2.23014e-12  8.34238e-16  1.0          2.25839e-12  2.25839e-12  8.34192e-16
 0.999998    1.21848e-6  1.21848e-6  4.44965e-9  1.0          2.23014e-12  2.23014e-12  8.34238e-16  1.0          2.25839e-12  2.25839e-12  8.34192e-16
 0.999998    1.21848e-6  1.21848e-6  4.44965e-9  1.0          2.23014e-12  2.23014e-12  8.34238e-16  1.0          2.25839e-12  2.25839e-12  8.34192e-16
 2.52546e-5  0.999973    2.03526e-7  1.21937e-6  0.000287185  0.999713     1.16563e-7   6.15819e-9   0.000315571  0.999684     1.1656e-7    6.76649e-9
 1.21848e-6  4.44965e-9  0.999998    1.21848e-6  2.23014e-12  8.34238e-16  1.0          2.23014e-12  2.25839e-12  8.34192e-16  1.0          2.25839e-12
 1.21848e-6  0.999998    4.44965e-9  1.21848e-6  2.23014e-12  1.0          8.34238e-16  2.23014e-12  2.25839e-12  1.0          8.34192e-16  2.25839e-12
 4.44965e-9  1.21848e-6  1.21848e-6  0.999998    8.34238e-16  2.23014e-12  2.23014e-12  1.0          8.34192e-16  2.25839e-12  2.25839e-12  1.0
 0.999998    1.21848e-6  1.21848e-6  4.44965e-9  1.0          2.23014e-12  2.23014e-12  8.34238e-16  1.0          2.25839e-12  2.25839e-12  8.34192e-16
 4.44965e-9  1.21848e-6  1.21848e-6  0.999998    8.34238e-16  2.23014e-12  2.23014e-12  1.0          8.34192e-16  2.25839e-12  2.25839e-12  1.0
 1.21848e-6  4.44965e-9  0.999998    1.21848e-6  2.23014e-12  8.34238e-16  1.0          2.23014e-12  2.25839e-12  8.34192e-16  1.0          2.25839e-12
 4.44965e-9  1.21848e-6  1.21848e-6  0.999998    8.34238e-16  2.23014e-12  2.23014e-12  1.0          8.34192e-16  2.25839e-12  2.25839e-12  1.0
 1.21848e-6  0.999998    4.44965e-9  1.21848e-6  2.23014e-12  1.0          8.34238e-16  2.23014e-12  2.25839e-12  1.0          8.34192e-16  2.25839e-12
 1.21848e-6  0.999998    4.44965e-9  1.21848e-6  2.23014e-12  1.0          8.34238e-16  2.23014e-12  2.25839e-12  1.0          8.34192e-16  2.25839e-12
 0.999998    1.21848e-6  1.21848e-6  4.44965e-9  1.0          2.23014e-12  2.23014e-12  8.34238e-16  1.0          2.25839e-12  2.25839e-12  8.34192e-16
 4.44965e-9  1.21848e-6  1.21848e-6  0.999998    8.34238e-16  2.23014e-12  2.23014e-12  1.0          8.34192e-16  2.25839e-12  2.25839e-12  1.0
 4.44965e-9  1.21848e-6  1.21848e-6  0.999998    8.34238e-16  2.23014e-12  2.23014e-12  1.0          8.34192e-16  2.25839e-12  2.25839e-12  1.0
]


# Another constructed data set
gene_mat=[[1 2 3 4 5 6]/sum([1 2 3 4 5 6]); [1 0 0 0 0 0]/sum([1 0 0 0 0 0]); [1001 1002 1003 1004 1005 1006]/sum([1001 1002 1003 1004 1005 1006]); 
    [1 2 3 4 5 6.1]/sum([1 2 3 4 5 6.1])]


# ## Check the data sets for different functions:

# +
# testing calckinship


function calckinship(geno::Matrix{Float64})

    # get dimensions
    sz = size(geno)

    # assign to variables for convenience
    nr = sz[1]
    nc = sz[2]

    # if empty then there is nothing to do
    if(nr==0)
        error("Nothing to do here.")
    else
        # make matrix to hold distances
        d = zeros(nr,nr)
    end

    # assign diagonals to ones
    for i=1:nr
        d[i,i] = 1.0
    end

    ncomplete = nc
    # off-diagonal elements need to be calculated
    if(nr>=2)
        for i=1:(nr-1)
            for j=(i+1):nr
                p1 = geno[i,:]
                p2 = geno[j,:]
                d[i,j] = d[j,i] = sum( p1 .* p2 + (1 .- p1) .* (1 .- p2) ) / ncomplete

            end
        end
    end
    return d
end
# -

calckinship(gpr)

calckinship(gene_mat)

# +
# testing kinship_man
function kinship_man(genematrix::Matrix{Float64})
#    c0=findall(.!isna.(genematrix[1,:]));
#    c1=findall(.!isna.(genematrix[2,:]));
#    ckeep=intersect(c0,c1);
#                  for j=3:nrow
#                  ck=find(.!isna.(genematrix[j,:]));
#                  ckeep=intersect(ck,ckeep);
#                   end
#                         geneupdate = genematrix[:,ckeep];
                        col = axes(genematrix,2)
                    kin=zeros(col,col);
                        @views for c=col, r=c:length(col)
                            kin[r,c]= 1.0-mean(abs.(genematrix[:,c]-genematrix[:,r]))
                                    kin[c,r]=kin[r,c]
                        end

return kin

end

# -

kinship_man(gpr)

kinship_man(gene_mat)

# testing kinship_lin
function kinship_lin(mat,cross)
r=size(mat,1)/cross; n=size(mat,2)
     K=Symmetric(BLAS.syrk('U','T',1.0,mat))/r
   @views for j=1:n
        K[j,j]=1.0
          end
    return convert(Array{Float64,2},K)
end


kinship_lin(gpr,1)

kinship_lin(gene_mat,1)

# +
# testing kinship_ctr
function kinship_ctr(genmat::Array{Float64,2})
   p=size(genmat,1)
    cgene= genmat.-mean(genmat,dims=2)
    # K=(1-ρ)*transpose(cgene)*cgene/n+ρ*Matrix(1.0I,n,n)
     #K=cgene'*cgene/p
     K=Symmetric(BLAS.syrk('U','T',1.0,cgene))/p
    return convert(Array{Float64,2},K)

end
# -

kinship_ctr(gpr)

kinship_ctr(gene_mat)

# testing kinship_std
function kinship_std(genmat::Array{Float64,2})
    p=size(genmat,1)
    sgene=(genmat.-mean(genmat,dims=2))./std(genmat,dims=2)
   #(1-ρ)*transpose(sgene)*sgene/n+ρ*Matrix(1.0I,n,n)
     K=Symmetric(BLAS.syrk('U','T',1.0,sgene))/p

    return convert(Array{Float64,2},K)
end

kinship_std(gpr)

kinship_std(gene_mat)

# +
# testing kinship_gs
function kinship_gs(climate::Array{Float64,2},ρ::Float64)
 env=axes(climate,2);
 Kc=zeros(env,env);

    @views for c=env, r=c:length(env)
            Kc[r,c]=exp(-mean(abs.(climate[:,c]-climate[:,r]).^2)/ρ)
            Kc[c,r]=Kc[r,c]
    end

    return Kc

end
# -

kinship_gs(gpr,.9)

kinship_gs(gene_mat,.9)

# +
# testing kinship_4way
function kinship_4way(genmat::Array{Float64,2})
    nmar, nind=axes(genmat)
    kmat=zeros(nind,nind);
    dist=zeros(nmar);

    for i=nind
        for j=i:length(nind)
            for k=nmar
                if (genmat[k,i]==genmat[k,j])
                    dist[k]= 0.0
                elseif (genmat[k,i]+genmat[k,j]==5)
                    dist[k]=2.0
                else
                    dist[k]=1.0
                end

            end
        kmat[j,i]=1.0-0.5*mean(dist)
        kmat[i,j]=kmat[j,i]
        end
    end

return kmat

end

# -

kinship_4way(gpr)

kinship_4way(gene_mat)

# ***All proporties of kinship are satisfied.***


