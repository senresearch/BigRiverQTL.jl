"""
parse_json(file::String)

Parses a JSON file to a dictionary with keys containing the names of CSV file for genomics data

# Argument

- `file` : A string containing the name of the JSON file

# Output

Returns a dictionary


"""

function parse_json(file::String)
	indict = open(file, "r") do f
		JSON.parse(f)
	end
	return indict
end





"""
read_data(filename::String)

Writes a CSV file to data frame excluding the comments lines.

# Argument

- `filename` : A string containing the name of the CSVfile

# Output

Returns a data frame.


"""
function read_data(filename)
	# read the file into lines
	lines = readlines(filename)
	# which lines have # as first character
	firstpound = (x->match(r"^#",x)).( lines ) .!= nothing
	# last line of comment
	startdata = findfirst(firstpound.==0)
    # Check if the comment lines can be run directly from "CSV.read". 
	return CSV.read(filename, DataFrame; header = startdata)#comment="#")
end

