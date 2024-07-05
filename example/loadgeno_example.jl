using BigRiverQTL
using BigRiverQTL.DataFrames
using BigRiverQTL.CSV

data_dir = joinpath(@__DIR__, "../data/BXD/")

file = joinpath(data_dir, "bxd.json")

bg_ctrl_dict = BigRiverQTL.parse_json(file)

keys(bg_ctrl_dict)

# get(bg_ctrl_dict, :geno)

genofile = joinpath(data_dir, bg_ctrl_dict["geno"])

function read_geno(filename)
        # read the file into lines
        lines = readlines(filename)
        # which lines have # as first character
        firstpound = (x->match(r"^#",x)).( lines ) .!= nothing
        # last line of comment
        endcomment = findfirst(firstpound.==0)
        
         return endcomment#CSV.read(filename, DataFrame; skipto = endcomment)

end



ecmt = read_geno(genofile)
CSV.read(genofile, DataFrame; skipto = ecmt)


str1 = "DXMit223,B,B,B,D,B,B,B,B,D,B,D,B,B,D,B,B,D,D,D,D,D,D,B,B,B,B,D,B,B,B,B,D,D,B,B,D,D,B,B,D,B,B,B,B,B,D,B,B,B,B,B,B,B,D,B,D,B,B,B,B,B,B,D,B,D,D,B,B,B,D,D,B,D,B,D,D,D,D,B,B,D,D,B,B,B,D,D,B,D,D,B,D,B,H,B,B,D,B,B,B,B,D,D,B,B,B,B,H,B,B,B,B,B,B,B,B,D,D,B,B,B,H,B,H,B,B,H,B,B,H,B,H,D,D,D,B,D,B,B,B,B,H,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,H,B,B,D,B,D,D,D,H,D,D,D,D,D,B,H,D,B,B,B,B,B,B,D,B,H,H,B,D,D,B,H,H,B"
str2 = "rs31638776,B,B,B,D,B,B,B,B,D,B,D,B,B,D,B,B,D,D,D,D,D,D,B,B,B,B,D,B,B,B,B,D,D,B,B,D,D,B,B,D,B,B,B,B,B,D,B,B,B,B,B,B,B,D,B,D,B,B,D,B,D,B,D,B,D,D,B,B,B,D,D,B,D,B,D,D,D,D,B,B,D,D,B,B,B,D,D,B,D,D,B,D,B,H,B,B,D,B,B,B,B,D,D,B,B,B,B,B,B,B,B,B,B,B,B,B,D,D,B,B,B,H,B,H,B,B,H,B,B,H,B,H,D,D,D,B,D,B,B,B,B,H,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,H,B,B,D,B,D,D,D,H,D,D,D,D,D,B,H,D,B,B,B,B,B,B,D,B,H,H,B,D,D,B,H,H,B"
str3 = "rs31639754,B,B,B,D,B,B,B,B,D,B,D,B,B,D,B,B,D,D,D,D,D,D,B,B,B,B,D,B,B,B,B,D,D,B,B,D,D,B,B,D,B,B,B,B,B,D,B,B,B,B,B,B,B,D,B,D,B,B,D,B,D,B,D,B,D,D,B,B,B,D,D,B,D,B,D,D,D,D,B,B,D,D,B,B,B,D,D,B,D,D,B,D,B,H,B,B,D,B,B,B,B,D,D,B,B,B,B,B,B,B,B,B,B,B,B,B,D,D,B,B,B,H,B,H,B,B,H,B,B,H,B,H,D,D,D,B,D,B,B,B,B,H,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,H,B,B,D,B,D,D,D,H,D,D,D,D,D,B,H,D,B,B,B,B,B,B,D,B,H,H,B,D,D,B,H,H,B"
str4 = "rs3693969,B,B,B,D,B,B,B,B,D,B,D,B,B,D,B,B,D,D,D,D,D,D,B,B,B,B,D,B,B,B,B,D,D,B,B,D,D,B,B,D,B,B,B,B,B,D,B,B,B,B,B,B,B,D,B,D,B,B,B,B,B,B,D,B,D,D,B,B,B,D,D,B,D,B,D,B,D,D,B,B,D,D,B,B,B,D,D,B,D,D,B,D,B,H,B,B,D,B,B,B,B,D,D,B,B,B,B,H,B,B,B,B,B,B,B,B,D,D,B,B,B,H,B,H,B,B,H,B,B,H,B,H,D,D,D,B,D,B,B,B,B,H,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,H,B,B,D,B,D,D,D,H,D,D,D,D,D,B,H,D,B,B,B,B,B,B,D,B,H,H,B,D,D,B,H,H,B"
split(str4, ",") |> length