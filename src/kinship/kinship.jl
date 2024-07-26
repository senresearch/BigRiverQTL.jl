#=
  `kinship.jl` contains functions for computing the genetic relatedness matrix (or kinship).
=#


"""
    calckinship(geno::Matrix{Float64})

Calculate kinship from the genotype probability array.

# Arguments

- `geno` is the genotype probability matrix for n individuals and p markers. 

# Output

Returns a n x n symmetric matrix containing 1's on the diagonal.

___

calckinship(geno::Matrix{Union{Missing,Float64}})

Calculate kinship from the genotype probability array.

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






