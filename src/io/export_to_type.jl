
"""
get_gmap(file::String)

Creates a `Gmap` type/struct from gmap CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the gmap CSV file.

# Output

Returns a `Gmap` type/struct.


"""
function get_gmap(filename::String)
    # load file   
    gdf = groupby(read_data(filename) , :chr)

    # make type
    
    # chromosome
    chr = [group.chr[1] for group in gdf];
    # marker
    marker = [String.(group.marker) for group in gdf]
    # relative position
    pos = [group.pos for group in gdf]

    return Gmap(chr, marker, pos)
end





"""
get_geno(file::String)

Creates a `Geno` type/struct from gmap CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the geno CSV file.

# Output

Returns a `Geno` type/struct.


"""
function get_geno(gmapfile::String,genofile::String)
    # load file   
    gmap=get_gmap(gmapfile);
    df_gmap=read_data(gmapfile)
    df_geno=read_data(genofile)

    # make type
    
    # samples
    samples=names(df_geno)[2:end]
    # chromosomes
    chr=gmap.chr
    # markers
    marker=gmap.marker
    # values
    val=[Matrix{Int}(undef,length(samples),length(marker[i])) for i in 1:length(chr)]
    for i in 1:length(chr)
        for j in 1:length(marker[i])
            uni=unique(Vector(df_geno[findfirst(x -> x==marker[i][j],marker[i]),2:end]))
           val[i][:,j]= recode(Vector(df_geno[findfirst(x -> x==marker[i][j],marker[i]),2:end]), uni[1]=>1, uni[2]=>2, uni[3]=>0)
        end
        
    end

    return Geno(samples,chr, marker, val)
end





"""
get_pmap(file::String)

Creates a `Pmap` type/struct from Pmap CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the pmap CSV file.

# Output

Returns a `Pmap` type/struct.


"""
function get_pmap(filename::String)
    # load file   
    gdf = groupby(read_data(filename) , :chr)

    # make type
    
    # chromosome
    chr = [group.chr[1] for group in gdf];
    # marker
    marker = [String.(group.marker) for group in gdf]
    # relative position
    pos = [group.pos for group in gdf]
    # unit
    unit= "mm10 Mbp"

    return Pmap(chr, marker, pos, unit)
end




"""
get_pheno(file::String)

Creates a `Pheno` type/struct from Pheno CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the pheno CSV file.

# Output

Returns a `Pheno` type/struct.


"""
function get_pheno(filename::String)
    # load file   
    
    df_pheno=read_data(filename)

    # make type
    
    # samples
    samples=df_pheno[:,1]

    # traits' names
    traits = names(df_pheno)[2:end]

    # value of the traits
    df_pheno=Matrix(df_pheno)
    val=tryparse.(Float64,df_pheno[:,2:end])

    return Pheno(samples,traits, val)
end




"""
get_phenocovar(file::String)

Creates a `Phenocovar` type/struct from phenocovar CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the phenocovar CSV file.

# Output

Returns a `Phenocovar` type/struct.


"""
function get_phenocovar(filename::String)
    # load file   
    
    df_phenocovar=read_data(filename)

    # make type
    
   # traits' names
   traits =string.(df_phenocovar[:,1])


   # description
   description=df_phenocovar[:,2]

    return Phenocov(traits, description)
end








"""
get_crossinfo(file::String)

Creates a `Crossinfo` type/struct from crossinfo CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the crossinfo CSV file.

# Output

Returns a `Crossinfo` type/struct.


"""
function get_crossinfo(filename::String)
    # load file   
    
    df_crossdirection=read_data(filename)

    # make type
    
   # traits' names
   samples =df_crossdirection[:,1]


   # description
   direction=df_crossdirection[:,2]

    return CrossInfo(samples, direction)
end










"""
get_isxchar(file::String)

Creates a `IsXChar` type/struct from gmap CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the gmap CSV file.

# Output

Returns a `IsXChar` type/struct.


"""
function get_isxchar(filename::String)
    # load file   
    
    gmap=get_gmap(filename)

    # make type

    # chromosome
    chr = gmap.chr

    # isXchar

    isxchar = zeros(Bool, length(chr))
    isxchar[findfirst(x->x=="X",chr)]=1

    return IsXChar(chr, isxchar)
end







"""
get_isfemale(file::String)

Creates a `IsFemale` type/struct from cross_info CSV file.

# Argument

- `filename` : A string containing the name(with directory) of the cross_info CSV file.

# Output

Returns a `IsFemale` type/struct.


"""
function get_isfemale(filename::String)
    # load file   
    
    jsondict=BigRiverQTL.parse_json(file)
    crossinfofile = joinpath(data_dir, jsondict["cross_info"]["file"])

    # make type

    # samples
    samples=get_crossinfo(crossinfofile).samples

    # isfemale

    isfemale = ones(Bool, length(samples))
    if(in("sex",keys(jsondict)))
        isfemale[findall(x->x=="Male",jsondict["sex"])]=0
    end

    return IsFemale(samples, isfemale)
end









"""
get_crosstype(file::String)

Creates a `CrossType` type/struct from  control file in json format.

# Argument

- `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `CrossType` type/struct.


"""
function get_crosstype(filename::String)
    # load file   
    
    jsondict=BigRiverQTL.parse_json(filename)

    # make type

    #crosstype
    if(in("crosstype",keys(jsondict)))
        crosstype=jsondict["crosstype"]

    else
        throw("Error: CrossType not found")
    end

    

    return CrossType(crosstype)
end









"""
get_alleles(file::String)

Creates a `Alleles` type/struct from  control file in json format.

# Argument

- `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `Alleles` type/struct.


"""
function get_alleles(filename::String)
    # load file   
    
    jsondict=BigRiverQTL.parse_json(filename)

    # make type

    #alleles
    alleles=jsondict["alleles"]

    

    return Alleles(alleles)
end










"""
get_genotype(file::String)

Creates a `GenoType` type/struct from  control file in json format.

# Argument

- `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `GenoType` type/struct.


"""
function get_genotype(filename::String)
    # load file   
    
    jsondict=BigRiverQTL.parse_json(filename)

    # make type

    #genotype
    if(in("genotypes",keys(bg_ctrl_dict)))
        label=jsondict["genotypes"]

    else
       label=Dict{String,Int}("A"=>1, "H"=>2, "B"=>3, "D"=>4, "C"=>5)
    end

    

    return GenoType(label)
end











"""
get_genotranspose(file::String)

Creates a `GenoTranspose` type/struct from  control file in json format.

# Argument

- `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `GenoTranspose` type/struct.


"""
function get_genotranspose(filename::String)
    # load file   
    
    jsondict=BigRiverQTL.parse_json(filename)

    # make type

    #genotype
    if(in("geno_transposed",keys(bg_ctrl_dict)))
        val=jsondict["geno_transposed"]

    else
       val=FALSE
    end

    

    return GenoTranspose(val)
end








"""
get_bigriverqtldata(file::String)

Creates a `BigRiverQTLData` type/struct from  control file in json format.

# Argument

- `filename` : A string containing the name(with directory) of the control file in json format.

# Output

Returns a `BigRiverQTLData` type/struct.


"""
function get_bigriverqtldata(filename::String)
    # load file   
    
    jsondict=BigRiverQTL.parse_json(filename)

    # make type

    #gmap
    if(in("gmap",keys(jsondict)))
        gmapfile = joinpath(data_dir, jsondict["gmap"])

    else
        throw("Error: gmap not found in control file")
    end
    
    gmap=get_gmap(gmapfile)

    #geno
    if(in("geno",keys(jsondict)))
        genofile = joinpath(data_dir, jsondict["geno"])

    else
        throw("Error: geno not found in control file")
    end
    
    geno=get_geno(gmapfile,genofile)

    #pmap
    if(in("pmap",keys(jsondict)))
        pmapfile = joinpath(data_dir, jsondict["pmap"])

    else
        throw("Error: pmap not found in control file")
    end
    
    pmap=get_pmap(pmapfile)

    #pheno
    if(in("pheno",keys(jsondict)))
        phenofile = joinpath(data_dir, jsondict["pheno"])

    else
        throw("Error: pheno not found in control file")
    end
    
    pheno=get_pheno(phenofile)

    #phenocov
    if(in("phenocovar",keys(jsondict)))
        phenocovfile = joinpath(data_dir, jsondict["phenocovar"])

    else
        throw("Error: phenocovar not found in control file")
    end
    
    phenocov=get_phenocovar(phenocovfile)

    #isXchar
    isXchar=get_isxchar(gmapfile)
    
    #isfemale
    isfemale=get_isfemale(filename)

    #crosstype
    crosstype=get_crosstype(filename)

    #crossinfo
    crossinfofile = joinpath(data_dir, jsondict["cross_info"]["file"])
    crossinfo=get_crossinfo(crossinfofile)
   
    #alleles
    alleles=get_alleles(filename)



    #genotype
    genotype=get_genotype(filename)



    #genotranspose
    genotranspose=get_genotranspose(filename)


    

    return BigRiverQTLData( gmap,
    geno,
    pmap,
    pheno,
    phenocov,
    isXchar,
    isfemale,
    crosstype,
    crossinfo,
    alleles,
    genotype,
    genotranspose)
end




