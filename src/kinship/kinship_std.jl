"""

     kinship_std(genmat::Array{Float64,2})


Calculates a kinship by a standardized (or normalized) genotype matrix (linear kernel), i.e. genotypes subtracted by marker mean and divided by marker standard deviation.
It can also do with climatic information data. See [`kinship_gs`](@ref).

# Argument

- `genmat` : A matrix of genotype data (0,1,2). size(genmat)=(n,p) for `n` individuals x `p` markers

# Output

Returns a n x n symmetric matrix.
See also [`kinship_ctr`](@ref).

"""
function kinship_std(genmat::Array{Float64,2})
    p=size(genmat,2)
    sgene=(genmat.-mean(genmat,dims=2))./std(genmat,dims=2)
   #(1-ρ)*transpose(sgene)*sgene/n+ρ*Matrix(1.0I,n,n)
     K=(sgene*sgene') ./ (p-1)

    return convert(Array{Float64,2},K)
end


