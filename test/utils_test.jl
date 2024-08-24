###############
# Get dataset #
###############

# get path
data_dir = joinpath(@__DIR__, "data/BXD/");
file = joinpath(data_dir, "bxd.json");

# Transforming data to a optimised and accessible data type
data = get_geneticstudydata(file);


########################
# Test encode_genotype #
########################

# Test Case 1: All genotypes are known
@testset "Encode Genotypes - Known Genotypes" begin
    geno_dict = Dict{String, Any}("AA" => 1, "AB" => 2, "BB" => 3)
    geno_val = ["AA" "AB" "BB"; "AB" "AA" "BB"]
    expected_output = [1 2 3; 2 1 3]
    @test encode_genotype(geno_dict, geno_val) == expected_output
end

# Test Case 2: Some genotypes are unknown
@testset "Encode Genotypes - Unknown Genotypes" begin
    geno_dict = Dict{String, Any}("AA" => 1, "AB" => 2, "BB" => 3)
    geno_val = ["AA" "AB" "BC"; "AB" "AA" "BD"]
    expected_output = [1 2 missing; 2 1 missing]
    encoded_output = encode_genotype(geno_dict, geno_val)
    @test (encoded_output[1:2,1:2] == expected_output[1:2,1:2]) && (ismissing(encoded_output[2,3]))
end

# Test Case 3: All genotypes are unknown
@testset "Encode Genotypes - All Unknown" begin
    geno_dict = Dict{String, Any}("AA" => 1, "AB" => 2)
    geno_val = ["AC" "AD"; "AE" "AF"]
    encoded_output = encode_genotype(geno_dict, geno_val)
    @test sum(ismissing.(encoded_output)) == 4
end


################
# Test missing #
################

test_geno = Geno(
    ["sample1", "sample2", "sample3"],
    ["1", "2"],
    [["marker1", "marker2", "marker3"], 
     ["marker4", "marker5", "marker6", "marker7"]],
    [[1 2 1; missing missing 2; 1 2 2],
     [1 2 2 1; missing missing 2 1; 2 1 2 missing]],
     CrossType("risib"),
     Alleles(["B", "D"]),
     GenoType(Dict{String, Any}("B" => 1, "D" => 2)),
     GenoTranspose(true)
);
test_geno_cleaned = get_geno_completecases(test_geno);

# Test the `get_geno_completecases` function
@testset "Testing get_geno_completecases function" begin
    @test typeof(test_geno_cleaned) == typeof(test_geno)  # Check type consistency
    @test sum(ismissing.(reduce(hcat, test_geno_cleaned.val))) == 0  # Ensure no missing values
    @test length(test_geno_cleaned.sample_id) == 2 # Check if some rows might have been removed
    @test length(test_geno_cleaned.marker_name[2]) == 3 # Check if some columns might have been removed
end

# Test the `summary_missing` function
@testset "Testing summary_missing function" begin
    tbl_missing = summary_missing(test_geno, issorted = true);

    @test tbl_missing[1].percentage[1] == 57.14 # Check missing percentage row wise
    @test tbl_missing.missing_per_col.percentage[1] == 33.33 # Check missing percentage column wise
end


# Test the `get_data_completecases` function
@testset "Testing get_data_completecases function" begin
    # Remove the  missing data
    data_no_missing = get_data_completecases(data);
    tbl_missing = summary_missing(data_no_missing.geno, issorted = true);

    @test sum(tbl_missing[1].percentage) == 0 # Check no missing row wise
    @test sum(tbl_missing[2].percentage) == 0 # Check no missing column wise
    @test data_no_missing.geno.sample_id == data_no_missing.pheno.sample_id # Ensure samples align  with pheno
    @test data_no_missing.geno.marker_name == data_no_missing.gmap.marker_name # Ensure loci align  with gmap
end

###############
# Test subset #
###############

# Test the `select_marker` function
@testset "Testing select_marker function" begin
    geno_subset_1 = select_marker(test_geno, ["marker2", "marker5"]);
    geno_subset_2 = select_marker(test_geno, Not(["marker3"]));
    geno_subset_3 = select_marker(test_geno, ["marker2"]);

    @test geno_subset_1.marker_name[1] == ["marker2"] # Check sample selection
    @test size(geno_subset_1.val[1], 2) == 1 # Check size of the geno matrix
    @test geno_subset_2.marker_name[1] == ["marker1", "marker2"] # Check sample selection
    @test size(geno_subset_2.val[1], 2) == 2 # Check size of the geno matrix
    @test geno_subset_3.marker_name[1] == ["marker2"] # Check sample selection
    @test size(geno_subset_3.val[1], 2) == 1 # Check size of the geno matrix
    @test length(geno_subset_3.chr) == 1 # Check size of the geno matrix
end


# Test the `select_sample` function
@testset "Testing select_sample function" begin
    geno_subset_1 = select_sample(test_geno, ["sample1"]);
    geno_subset_2 = select_sample(test_geno, Not(["sample3"]));

    @test geno_subset_1.sample_id == ["sample1"] # Check sample selection
    @test size(geno_subset_1.val[1], 1) == 1 # Check size of the geno matrix
    @test geno_subset_2.sample_id == ["sample1", "sample2"] # Check sample selection
    @test size(geno_subset_2.val[1], 1) == 2 # Check size of the geno matrix
end



###################
# Test converting #
###################

test_gmap = Gmap(
    ["1", "2"],
    [["marker1", "marker2", "marker3"], 
     ["marker4", "marker5", "marker6", "marker7"]],
     [[1.58, 2.75, 3.59], 
     [0.58, 3.67, 8.87, 12.03]],
     "cM"
);

# Test the `gmap2df` function
@testset "Testing gmap2df function" begin
    df_g = gmap2df(test_gmap);
    
    @test isa(df_g, DataFrame) # Check type of output
    @test size(df_g) == (7, 3) # Check size of the gmap dataframe
end

test_pmap = Pmap(
    ["1", "2"],
    [["marker1", "marker2", "marker3"], 
     ["marker4", "marker5", "marker6", "marker7"]],
     [[1.58, 2.75, 3.59], 
     [0.58, 3.67, 8.87, 12.03]],
     "cM"
);


# Test the `pmap2df` function
@testset "Testing pmap2df function" begin
    df_p = pmap2df(test_pmap);
    
    @test isa(df_p, DataFrame) # Check type of output
    @test size(df_p) == (7, 3) # Check size of the gmap dataframe
end
