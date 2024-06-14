#=
  `kinship.jl` contains functions for computing genetic relatedness matrix (or kinship).
=#


"""
    calckinship(geno::Matrix{Float64})

Calculate kinship from genotype probability array.

# Arguments

- `geno` is the genotype probability matrix for n individuals and p markers. 

# Output

Returns a n x n symmetric matrix containing 1's on the diagonal.

___

calckinship(geno::Matrix{Union{Missing,Float64}})

Calculate kinship from genotype probability array.

# Arguments

- `geno` is the genotype probability matrix, n individuals and p markers, which 
    contains `missing` values.

# Output

Returns a n x n symmetric matrix containing 1's on the diagonal.

"""
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


function calckinship(geno::Matrix{Union{Missing,Float64}})

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

    iscomplete = Array{Bool,1}(undef,nc)
    ncomplete::Int64 = 0
    # off-diagonal elements need to be calculated
    if(nr>=2)
        for i=1:(nr-1)
            for j=(i+1):nr
                iscomplete = .!( ismissing.(geno[i,:]) .& ismissing.(geno[j,:]) )
                ncomplete = sum(iscomplete)
                p1 = geno[i,iscomplete]
                p2 = geno[j,iscomplete]
                d[i,j] = d[j,i] = sum( p1 .* p2
                                       + (1-p1) .* (1-p2) ) / ncomplete

            end
        end
    end
    return d
end



"""

      kinshipMan(genematrix::Matrix{Float64})

Calculates a kinship matrix using a manhattan distance. Missing values need to be either omitted or imputed.
This function is for recombinant inbred line (RIL) (AA/BB), not for 4-way cross genotype data.  See [`kinship4way`](@ref).

# Argument

- `genematrix` : A matrix of genotypes, i.e. 0,1 (or 1,2).  size(genematrix)= (p,n) for `p` genetic markers x `n` individuals(or lines).


# Output

Returns a n x n symmetric matrix containing 1's on the diagonal.

"""
function kinshipMan(genematrix::Matrix{Float64})
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
