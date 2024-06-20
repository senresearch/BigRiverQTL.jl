"""

     kinship_ctr(genmat::Array{Float64,2})

Calculates a kinship by a centered genotype matrix (linear kernel), i.e. genotypes subtracted by marker mean.

# Argument

- `genmat` : A matrix of genotype data (0,1,2). size(genmat)=(p,n) for `p` markers x `n` individuals

# Output

Returns a n x n symmetric matrix.
See also [`kinship_std`](@ref).

"""
function kinship_ctr(genmat::Array{Float64,2})
   p=size(genmat,1)
    cgene= genmat.-mean(genmat,dims=2)
    # K=(1-ρ)*transpose(cgene)*cgene/n+ρ*Matrix(1.0I,n,n)
     #K=cgene'*cgene/p
     K=Symmetric(BLAS.syrk('U','T',1.0,cgene))/p
    return convert(Array{Float64,2},K)

end
