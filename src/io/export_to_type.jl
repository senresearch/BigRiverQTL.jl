
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
    df_pheno[findall(x->x=="NA",df_pheno)].="-9999"
    val=parse.(Float64,df_pheno[:,2:end])

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
   traits =df_phenocovar[:,1]


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
   samples =df_phenocovar[:,1]


   # description
   direction=df_phenocovar[:,2]

    return CrossInfo(samples, direction)
end