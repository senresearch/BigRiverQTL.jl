module BigRiverQTL
    using BulkLMM
    using DataFrames

    include("loco_helpers.jl")
    export get_loco_geno, get_loco_geno_info, calcLocoKinship

    include("loco_scan.jl")
    export loco_scan

end
