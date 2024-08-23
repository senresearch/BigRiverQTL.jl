"""

     shrinkg(f,nb::Int64,geno)

Estimates a full-rank positive definite kinship matrix by shrinkage intensity estimation (bootstrap).  It can only be used with [`kinship_man`](@ref), [`kinship_4way`](@ref).
This function runs faster by CPU parallelization.  Add workers/processes using `addprocs()` function before running for speedup.

# Arguments

* `f `: A function of computing a kinship. It can only be used with [`kinship_man`](@ref), [`kinship_4way`](@ref).
* `nb` : An integer indicating the number of bootstraps. It does not have to be a large number.
* `geno` : A matrix of genotypes. See [`kinship_man`](@ref), [`kinship_4way`](@ref) for dimension.

# Example

```julia
julia> using flxQTL
julia> addprocs(8)
julia> K = shinkage(kinshipMan,20,myGeno)
```

# Output
_
Returns a full-rank symmetric positive definite matrix.

"""
function shrinkg(f,nb::Int64,geno)
    # generate a kinship matrix
    K0=f(geno)
    n=size(K0,1);
    # a target matrix
#             Id=Matrix(1.0I,n,n)
    #compute an optimal regularization parameter λ_hat
    denom=norm(I-K0)^2
    #np=nprocs();
    GG=[];

    idx=pmap(sample,[1:n for i=1:nb],[n for i=1:nb])
    Genematrix=pmap(f,[geno[:,idx[i]] for i=1:nb])
    GG=[GG;Genematrix]

    Ks=zeros(n,n,nb)
   @views for t=1:nb
        Ks[:,:,t]=GG[t]
    end
    kinVar=var(Ks;dims=3);
λ_hat=sum(kinVar)/denom;
K=λ_hat*I+(1-λ_hat)*K0;

# println("λ_hat is"," ",λ_hat,".")
return K
end

# function shrinkg(f,nb,cross,geno)
#     # generate a kinship matrix
#     K0=f(geno,cross)
#     n=size(K0,1);
#     # a target matrix
# #             Id=Matrix(1.0I,n,n)
#     #compute an optimal regularization parameter λ_hat
#     denom=norm(I-K0)^2
#    #np=nprocs(); itr=fld(nb,np)  #floor(nb/np);
#     GG=[];

#     idx=pmap(sample,[1:n for i=1:nb],[n for i=1:nb])
#     Genematrix=pmap(f,[geno[:,idx[i]] for i=1:nb],[cross for i=1:nb])
#     GG=[GG;Genematrix]

#     Ks=zeros(n,n,nb)
#    @views for t=1:nb
#         Ks[:,:,t]=GG[t]
#     end
#     kinVar=var(Ks;dims=3);
# λ_hat=sum(kinVar)/denom;
# K=λ_hat*I+(1-λ_hat)*K0;
# # println("λ_hat is"," ",λ_hat,".")
# return K
# end
