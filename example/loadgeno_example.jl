using BigRiverQTL
using BigRiverQTL.DataFrames
using BigRiverQTL.CSV

data_dir = joinpath(@__DIR__, "../data/BXD/")

file = joinpath(data_dir, "bxd.json")

bg_ctrl_dict = BigRiverQTL.parse_json(file)

keys(bg_ctrl_dict)

# get(bg_ctrl_dict, :geno)

genofile = joinpath(data_dir, bg_ctrl_dict["geno"])

function read_data(filename)
        # read the file into lines
        lines = readlines(filename)
        # which lines have # as first character
        firstpound = (x->match(r"^#",x)).( lines ) .!= nothing
        # last line of comment
        startdata = findfirst(firstpound.==0)

        return CSV.read(filename, DataFrame; header=startdata)
end

gmapfile = joinpath(data_dir, bg_ctrl_dict["gmap"])


function get_gmap(filename)
    # load file   
    gdf = groupby(read_data(filename) , :chr)
    
    chr = [group.chr[1] for group in gdf];
    marker = [String.(group.marker) for group in gdf]
    pos = [group.pos for group in gdf]

    return Gmap(chr, marker, pos)
end

gmap = get_gmap(gmapfile)


df_test = read_data(gmapfile)
gdf_test = groupby(df_test , :chr)


marker = [String.(group.marker) for group in gdf]


gdf_test[20].pos[1]== gdf_test[20].pos[2];



genofile = joinpath(data_dir, bg_ctrl_dict["geno"])
df_gmap = read_gmap()
df_geno = read_geno(genofile)
