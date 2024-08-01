module BigRiverQTL
    using BulkLMM
    using DataFrames, JSON, CSV, CategoricalArrays
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
    include("struct/datastructure.jl")
    export Gmap, Alleles, CrossType, GenoType, GenoTranspose,Geno, Pmap, Pheno, Phenocov, IsFemale, IsXChar,  CrossInfo
    export GeneticStudyData

    ######
    # IO #
    ######
    include("io/io_utils.jl")
    include("io/export_to_type.jl")
    export get_geneticstudydata


    ##########################
    # using BigRiverQTLPlots #
    ##########################
    using Reexport
    @reexport import BigRiverQTLPlots: plot_QTL, plot_eQTL, plot_manhattan

    #########
    # Plots #
    #########
    include("plots/plots_utils.jl")
    export gmap2df pmap2df

    include("plots/plots_qtl.jl")
    export plot_QTL 

    include("plots/plots_manhattan.jl")
    export plot_manhattan


    include("plots/plots_eqtl.jl")
    export plot_eQTL 






    
    






    
    

end
