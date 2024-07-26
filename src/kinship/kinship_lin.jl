"""

    kinship_lin(mat,cross)

Calculates a kinship (or climatic relatedness, [`kinship_gs`](@ref)) matrix by linear kernel.

# Arguments

- `mat` : A matrix of genotype (or allele) probabilities usually extracted from [R/qtl](https://rqtl.org/tutorials/rqtltour.pdf),
        [R/qtl2](https://kbroman.org/qtl2/assets/vignettes/user_guide.html) or the counterpart packages. size(mat)= (p,n) for p genetic markers x n individuals.
- `cross` : A scalar indicating instances of alleles or genotypes in a genetic marker. ex. 1 for genotypes (labeled as 0,1,2), 2 for RIF, 4 for four-way cross, 8 for HS mouse (allele probabilities), etc.

# Output

Returns a n x n symmetric (positive definite) matrix containing 1's on the diagonal.

See also [`kinshipCtr`](@ref), [`kinshipStd`](@ref) for genetype data.


"""
function kinship_lin(mat,cross)
r=size(mat,1)/cross; n=size(mat,2)
     K=Symmetric(BLAS.syrk('U','T',1.0,mat))/r
   @views for j=1:n
        K[j,j]=1.0
          end
    return convert(Array{Float64,2},K)
end