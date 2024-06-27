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

# # Testing Kinship functions

using Pkg, Revise

pwd()

Pkg.activate("..")

Pkg.instantiate()

# ### Libraries

using BigRiverQTL

using Statistics
using Distributed
using LinearAlgebra
import StatsBase: sample
# include("Miscellanea.jl")

# ### Data

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


# Another constructed data set for functions like "kinship_ctr" and "kinship_std".
gene_mat2=rand([0.0, 1.0, 2.0], 5, 10)

# ## Check the data sets for different functions:

# ### We will check the following properties for each kinship matrix: 
# * Dimensions: $n \times n$ where $n$ is the number of subjects.
# * All elements in the diagonals are $1$.
# * All elements in each kinship matrix lie between $-1$ and $1$.
# * Symmetric
# * Positive Definite
#
# We create a function "kinship_test" to test all the mentioned properties. 

# +
# Checks the above mentioned properties
function kinship_test(kinmat::Matrix{Float64})
    a=zeros(4)
    
    #checking dimensions
    if(size(kinmat,1)==size(kinmat,2))
        a[1]=1
        else(a[1]==0)
        println("The matrix is not a square matrix.")
    end
    
    # checking if all elements between -1 and 1
    if(maximum(kinmat)<=1.0 && minimum(kinmat)>=-1.0)
        a[2]=1
        else(a[2]==0)
        println("Not all elements of the matrix lie between -1 and 1.")
    end

    # checking if all elements in the diagonal are 1
    if(diag(kinmat)==ones(size(kinmat,1)))
        a[3]=1
        else(a[3]==0)
        println("Not all elements of the diagonal  are 1.")
    end

    # checking if 'kinmat' symmetric and positive definite
    if(isposdef(kinmat)==1)
        a[4]=1
        else(a[4]==0)
        println("The matrix is not positive definite.")
    end
    
    if(a==ones(4))
        println("All mentioned properties of a kinship matrix are satisfied.")
        
    
    
    
    
    
    
        
    
        

    end
        
        
end
    
# -

# ## Testing "calckinship"

A=calckinship(gpr)

kinship_test(A)

A=calckinship(gene_mat)

kinship_test(A)

# ***Testing "kinship_man"***

A=kinship_man(gpr)

kinship_test(A)

A=kinship_man(gene_mat)

kinship_test(A)

# ***Testing "kinship_lin"***


A=kinship_lin(gpr,1)

kinship_test(A)

kinship_lin(gene_mat,1)

kinship_test(A)

# ***Testing "kinship_ctr"***

A=kinship_ctr(gene_mat2)

kinship_test(A)





# ***Testing "kinship_std"***

A=kinship_std(gene_mat2)

kinship_test(A)

# ***Testing "kinship_gs"***

A=kinship_gs(gpr,.9)

kinship_test(A)





eigvals(A)

A=kinship_gs(gene_mat,.9)

kinship_test(A)

# ***Testing "kinship_4way"***

A=kinship_4way(gpr)

kinship_test(A)

A=kinship_4way(gene_mat)

kinship_test(A)

genmat=gene_mat2.-1
p=size(gene_mat2,2)

cgene= genmat#.-mean(genmat,dims=2)
    # K=(1-ρ)*transpose(cgene)*cgene/n+ρ*Matrix(1.0I,n,n)
     #K=cgene'*cgene/p
     K=cgene'*cgene/p

eigvals(K)

det(K)

K|>x->eigvals(x)|>x->round.(x,digits=15)
