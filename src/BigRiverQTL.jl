module BigRiverQTL
    using BulkLMM
    using DataFrames

    ########
    # Loco #
    ########

    include("loco/loco_helpers.jl")
    export get_loco_geno, get_loco_geno_info, calcLocoKinship

    include("loco/loco_scan.jl")
    export loco_scan

    include("loco/loco_bulkscan.jl")
    export loco_bulkscan

    ###########
    # Kinship #
    ###########




end
