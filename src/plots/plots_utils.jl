function gmap2df(gmap::Gmap)
    vmarker=reduce(vcat,gInfo.marker_name);
    Mb=reduce(vcat,gInfo.pos);
    start=0
    vchr=repeat(["0"],sum(length.(gInfo.marker_name)))
    for i in 1: eachindex(gInfo.chr)
        l_i=length(gInfo.marker_name[i])
        vchr[start+1:l_i+start] .= gInfo.chr[i]
        start=start+l_i
    end


end