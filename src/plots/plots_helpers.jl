function gmap2df(gmap::Union{Gmap, Pmap})
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


function gmap2df(gmap::Gmap)
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


function pmap2df(pmap::Pmap)

    start=0
    Chr =repeat(["0"],sum(length.(pmap.marker_name)))
        
    for i in eachindex(pmap.chromosomes)
        l_i=length(pmap.marker_name[i])
        Chr[start+1:l_i+start] .= pmap.chr[i]
        start=start+l_i
    end
    
    df=DataFrame(
        Locus = reduce(vcat,pmap.marker_name),
        Chr = Chr,
        Pos= reduce(vcat,pmap.pos)
    )
    return df
end