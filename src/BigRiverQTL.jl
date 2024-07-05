module BigRiverQTL
    using BulkLMM
    using DataFrames, JSON, CSV
    using Statistics
    using Distributed
    using LinearAlgebra
    import StatsBase: sample
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

    include("kinship/kinship.jl")
    export calckinship

    include("kinship/kinship_4way.jl")
    export kinship_4way

    include("kinship/kinship_ctr.jl")
    export kinship_ctr

    include("kinship/kinship_gs.jl")
    export kinship_gs

    include("kinship/kinship_lin.jl")
    export kinship_lin

    include("kinship/kinship_man.jl")
    export kinship_man

    include("kinship/kinship_std.jl")
    export kinship_std

    include("kinship/shrinkg.jl")
    export shrinkg

    #############
    # Structure #
    #############

    include("io/parse_json.jl")

end
