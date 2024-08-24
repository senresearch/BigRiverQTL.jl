"""
    map2df(gmap::Union{Gmap, Pmap}) -> DataFrame

Convert a genetic or physical map object to a DataFrame.

# Arguments
* `gmap`: A genetic map (`Gmap`) or phenotype map (`Pmap`) object that contains mapping information. 
  The object is expected to have at least three fields:
  * `chr`: An array where each element is a chromosome identifier corresponding to markers.
  * `marker_name`: A nested array where each sub-array contains the names of markers for a corresponding chromosome.
  * `pos`: A nested array similar to `marker_name`, but containing the positions of markers.

# Returns
* `DataFrame`: A DataFrame with three columns:
  * `Locus`: A flat list of all marker names.
  * `Chr`: A list of chromosome identifiers corresponding to each marker.
  * `Pos`: A flat list of all marker positions.

"""
function map2df(gmap::Union{Gmap, Pmap})
    start=0
    Chr =repeat(["0"],sum(length.(gmap.marker_name)))
       
    for i in eachindex(gmap.chr)
        l_i=length(gmap.marker_name[i])
        Chr[start+1:l_i+start] .= gmap.chr[i]
        start=start+l_i
    end
    df=DataFrame(
        Locus = reduce(vcat,gmap.marker_name),
        Chr = Chr,
        Pos= reduce(vcat,gmap.pos)
    )
    return df
end


"""
    gmap2df(gmap::Gmap) -> DataFrame

Convert a genetic map object (`Gmap`) into a DataFrame.

# Arguments
- `gmap`: A genetic map object containing mapping information. 
  This object should conform to the expected structure, which includes:
  - `chr`: An array of chromosome identifiers.
  - `marker_name`: A nested array of marker names organized by chromosome.
  - `pos`: A nested array of marker positions, also organized by chromosome.

# Returns
- `DataFrame`: A DataFrame representing the genetic map data, with columns 
for marker names (`Locus`), chromosome identifiers (`Chr`), and 
marker positions (`Pos`).

"""
function gmap2df(gmap::Gmap)
    return map2df(gmap)
end


"""
    pmap2df(pmap::Pmap) -> DataFrame

Convert a phenotype map object (`Pmap`) into a DataFrame.

# Arguments
- `pmap`: A genetic map object containing mapping information. 
  This object should conform to the expected structure, which includes:
  - `chr`: An array of chromosome identifiers.
  - `marker_name`: A nested array of marker names organized by chromosome.
  - `pos`: A nested array of marker positions, also organized by chromosome.

# Returns
- `DataFrame`: A DataFrame representing the phenotypic map data, with columns 
for marker names (`Locus`), chromosome identifiers (`Chr`), and 
marker positions (`Pos`).

"""
function pmap2df(pmap::Pmap)
    return map2df(pmap)
end