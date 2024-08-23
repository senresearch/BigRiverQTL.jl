"""
plot_manhattan(vLOD::Vector{<: AbstractFloat}, dfgInfo::DataFrame;
		chrColname::String="Chr", mbColname::String="Mb",
		thresholds::Vector{<: AbstractFloat}=[], kwargs...)

Generates a scatter plot for QTL analysis.

## Arguments
* `vLOD` is the vector containing the maximum value of LOD score.
* `dfgInfo` is a data of type `Gmap` containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
* `chrColname` is the name of the column containing the chromosomes' names, default name is "Chr".
* `mbColname` is the name of the column containing the megabase DNA length, default name is "Mb". 
* `thresholds` is <: AbstractFloat number vector containing desired LOD score thresholds for plotting.

---

plot_manhattan(scanresult::NamedTuple, dfgInfo::DataFrame;
		chrColname::String="Chr", mbColname::String="Mb", 
		thresholds::Vector{<: AbstractFloat}=[], kwargs...)

Generates a scatter plot for QTL analysis.

## Arguments
* `scanresult` is NamedTuple object resulting from the `scan()` function in `BulkLMM.jl`.
* `dfgInfo` is a data of type `Gmap` containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
* `chrColname` is the name of the column containing the chromosomes' names, default name is "Chr".
* `mbColname` is the name of the column containing the megabase DNA length, default name is "Mb". 
* `significance` is <: AbstractFloat number vector containing significant levels to estimate LOD score thresholds.

If the `scanresult` does not contain a permutation matrix, the original maximum LOD scores will be plotted, and the values in 
the `significance` vector will be used as the threshold values for comparison.

"""
function plot_manhattan(vLOD::Vector{<: AbstractFloat}, dfgInfo::Gmap;
	kwargs...)

	return plot_manhattan(vLOD, gmap2df(dfgInfo); kwargs...)

end


function plot_manhattan(scanresult::NamedTuple, dfgInfo::Gmap;
	kwargs...)

	return plot_manhattan(scanresult, gmap2df(dfgInfo); kwargs...)

end





function plot_manhattan!(vLOD::Vector{<: AbstractFloat}, dfgInfo::Gmap;
	kwargs...)

	return plot_manhattan!(vLOD, gmap2df(dfgInfo); kwargs...)

end


function plot_manhattan!(scanresult::NamedTuple, dfgInfo::Gmap;
	kwargs...)

	return plot_manhattan!(scanresult, gmap2df(dfgInfo); kwargs...)

end







