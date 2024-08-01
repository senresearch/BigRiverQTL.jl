
"""
plot_eQTL(multiLODs::Array{Float64, 2}, dfpInfo::DataFrame, dfgInfo::DataFrame;
         threshold::Float64 = 5.0, kwargs...)

Generates a scatter plot for eQTL analysis.

## Arguments
- `multiLODs` is the matrix containing the LOD's values.
- `dfpInfo` is a data of type `Gmap` containing the phenotype information such as probeset, chromosomes names and Mb distance.
- `dfgInfo` is a data of type `Pmap` containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
- `threshold` is the LOD threshold value, default is `5.0``.

"""
function plot_eQTL(multiLODs::Array{Float64, 2}, dfpInfo::Pmap, dfgInfo::Gmap;
                   kwargs...)

                   return plot_eQTL(multiLODs,pmap2df(dfpInfo), gmap2df(dfgInfo); kwargs...)

end
