#The genotype data can be extracted from `cross object` using [r-qtl](https://rqtl.org/tutorials/) or [r-qtl2](https://kbroman.org/qtl2/assets/vignettes/user_guide.html).

"""

     kinship_4way(genmat::Array{Float64,2})

Computes a kinship for four-way cross data counting different alleles between two markers: ex. AB-AB=0; AB-AC=1; AB-CD=2,``\\dots``
Note: In [R/qtl](https://cran.r-project.org/web/packages/qtl/qtl.pdf), genotypes are labeled as 1=AC; 2=BC; 3=AD; 4=BD by the function, `read.cross`.


# Argument

* `genmat` : A matrix of genotypes for `four-way cross` ``(1,2, \\dots)``.
           size(genematrix)= (p,n), for `p` genetic markers x `n` individuals(or lines).

# Output

Returns a n x n symmetric matrix containing 1's on the diagonal.

"""
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


