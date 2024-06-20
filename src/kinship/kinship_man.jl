"""

      kinship_man(genematrix::Matrix{Float64})

Calculates a kinship matrix using a Manhattan distance. Missing values need to be either omitted or imputed.
This function is for recombinant inbred line (RIL) (AA/BB), not for 4-way cross-genotype data.  See [`kinship_4way`](@ref).

# Argument

- `genematrix` : A matrix of genotypes, i.e. 0,1 (or 1,2).  size(genematrix)= (p,n) for `p` genetic markers x `n` individuals(or lines).


# Output

Returns a n x n symmetric matrix containing 1's on the diagonal.

"""
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
