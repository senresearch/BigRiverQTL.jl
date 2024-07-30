function gmap2df(gmap::Gmap)
    Locus=reduce(vcat,gmap.marker_name);
    Mb=reduce(vcat,gmap.pos);
    start=0
    Chr=repeat(["0"],sum(length.(gmap.marker_name)))
    Cm=repeat([Missing],sum(length.(gmap.marker_name)))
    for i in  eachindex(gmap.chr)
        l_i=length(gmap.marker_name[i])
        Chr[start+1:l_i+start] .= gmap.chr[i]
        start=start+l_i
    end
    df=DataFrame(Locus=Locus,Chr=Chr,Mb=Mb,cM=Cm)
    return df

end