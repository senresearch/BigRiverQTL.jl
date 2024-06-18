"""

     kinshipStd(genmat::Array{Float64,2})


Calculates a kinship by a standardized (or normalized) genotype matrix (linear kernel), i.e. genotypes subtracted by marker mean and divided by marker standard deviation.
Can also do with climatic information data. See [`kinshipGs`](@ref).

# Argument

- `genmat` : A matrix of genotype data (0,1,2). size(genmat)=(p,n) for `p` markers x `n` individuals

# Output

Returns a n x n symmetric matrix.
See also [`kinshipCtr`](@ref).

"""
function kinship_std(genmat::Array{Float64,2})
    p=size(genmat,1)
    sgene=(genmat.-mean(genmat,dims=2))./std(genmat,dims=2)
   #(1-ρ)*transpose(sgene)*sgene/n+ρ*Matrix(1.0I,n,n)
     K=Symmetric(BLAS.syrk('U','T',1.0,sgene))/p

    return convert(Array{Float64,2},K)
end


