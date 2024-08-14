
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
