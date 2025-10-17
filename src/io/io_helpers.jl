"""
parse_json(file::String)

Parses a JSON file to a dictionary with keys containing the names of CSV file for genomics data

# Argument

* `file` : A string containing the name of the JSON file

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

* `filename` : A string containing the name of the CSVfile

# Output

Returns a data frame.
"""
function read_data(filename::String; kwargs...)
	# # read the file into lines
	# lines = readlines(filename)
	# # which lines have # as first character
	# firstpound = (x->match(r"^#",x)).( lines ) .!= nothing
	# # last line of comment
	# startdata = findfirst(firstpound.==0)
    # # Check if the comment lines can be run directly from "CSV.read". 
	return CSV.read(filename, DataFrame; kwargs...)#comment="#")
end



"""
    get_control_file(filename::String) -> (String, String)

Locate the control file and determine the directory where it resides.

# Arguments
* `filename::String`: A string representing either the path to a file or a directory.

# Returns
* `(String, String)`: A tuple containing the directory of the control file and the path 
to the control file. If the provided `filename` is a directory, the function will search 
for JSON files within and return the path to the first JSON file found.

# Description
The function first resolves the absolute path of the provided `filename`. If `filename` 
points to an existing file, the function identifies the directory containing this 
file. If `filename` specifies a directory, the function lists all entries in this 
directory, searches for the first JSON file, and updates `filename` to the path of this 
JSON file.
"""
function get_control_file(filename::String)

    filename = realpath(filename)
    if isfile(filename)
        data_dir = dirname(realpath(filename))
    elseif isdir(filename)
        data_dir = filename
        files_dir = readdir(data_dir) 
        filename = files_dir |> 
                    x ->  findall(occursin.(".json", filesdir)) |>
                    x -> filesdir[x] |>
                    x -> joinpath(data_dir, x) 
    end

	return data_dir, filename
end

"""
    encode_genotype(geno_dict::Dict{String, Any}, geno_val::AbstractArray) -> Matrix

Encode a matrix of genotype values using a dictionary of predefined mappings.

# Arguments
* `geno_dict::Dict{String, Any}`: A dictionary mapping genotype strings to integer codes.
* `geno_val::AbstractArray`: An array of genotype strings to be encoded.

# Returns
* `Matrix{Union{Missing, Int64}}` or `Matrix{Int64}`: A matrix of the same dimensions as 
`geno_val` where each genotype string is replaced by its corresponding integer code. 
Genotypes not found in `geno_dict` are encoded as `missing`.

# Description
The function first copies the provided `geno_dict` to preserve the original. It then identifies 
any genotypes in `geno_val` that are not already present in the dictionary. A warning is issued 
for these missing genotypes, and they are added to the dictionary with a mapping to `missing`.

Each element in the input `geno_val` is replaced in the `encoded_val` matrix by its corresponding 
integer code from `geno_dict`, handling both existing and newly added genotype mappings.

# Examples
```julia
geno_dict = Dict{String, Any}("AA" => 1, "AB" => 2, "BB" => 3)
geno_val = ["AA", "AB", "BB", "BC"]
encoded_matrix = encode_genotype(geno_dict, geno_val)
println(encoded_matrix)
```
"""
function encode_genotype(geno_dict::Dict{String, Any}, geno_val::AbstractArray)
    # check if geno_val is a Vector
    if size(geno_val, 2) == 1
        geno_val = reshape(geno_val, size(geno_val, 1), 1)
    end

	# need to copy dict to preserve original
	geno_dict = copy(geno_dict)
	
    # check uniqueness of genotype
	unique_genotype = unique(geno_val)
	dict_genotype = string.(keys(geno_dict))
	extra_genotype = setdiff(unique_genotype, dict_genotype)

    # initialize output
	if !isempty(extra_genotype)
		@warn "Presence of genotypes not encoded by dictionnary of the control file: " * join(extra_genotype, " & ") *
			  "\n Not specified genetotype will be encoded as missing"

		for i in eachindex(extra_genotype)
			setindex!(geno_dict, missing, extra_genotype[i])
		end
		encoded_val = Array{Union{Missing, Int64}}(undef, size(geno_val))
	else
		encoded_val = Array{Int64}(undef, size(geno_val))
	end

    # encoding
	for j in 1:size(geno_val, 2)
		for i in 1:size(geno_val, 1)
			encoded_val[i, j] = geno_dict[geno_val[i, j]]
		end
	end

	return encoded_val
end


"""
    check_key(control_dict::Dict, s::String) -> Any

Check if a specified key exists in the given dictionary and return its corresponding value.

# Arguments
* `control_dict::Dict`: A dictionary from which the value associated with a key is to be retrieved.
* `s::String`: The key for which the existence and value are checked within the dictionary.

# Returns
* `Any`: Returns the value associated with the key `s` in `control_dict` if it exists.

# Throws
- Throws an error if the key `s` is not found in `control_dict`.

# Description
This function checks if the given string `s` exists as a key in the provided dictionary 
`control_dict`. If the key exists, the function returns the value associated with this 
key. If the key does not exist, it throws an error with a message indicating that the 
key was not found in the control file.

"""
function check_key(control_dict::Union{Dict, JSON.Object}, s::String)
	# check geno file exists
	if (in(s, keys(control_dict)))
		val = control_dict[s]
	else
		@warn "$(s) file not found in control file"
		val = missing
		# throw("Error: $(s) not found in control file")
	end
	
	return val
end




