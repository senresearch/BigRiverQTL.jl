module BigRiverQTL
    using BulkLMM
    using DataFrames, JSON, CSV
    using Statistics
    using LinearAlgebra
    import StatsBase: sample

    using Reexport
    @reexport import BigRiverQTLPlots: plot_QTL, plot_eQTL, plot_manhattan
    @reexport import BulkLMM: scan


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
    export Gmap, Alleles, CrossType, GenoType, GenoTranspose, Geno, Pmap
    export Pheno, Phenocov, IsFemale, IsXChar, CrossInfo
    export GeneticStudyData

    ######
    # IO #
    ######
    include("io/io_utils.jl")
    export get_control_file, encode_genotype
    include("io/export_to_type.jl")
    export get_geneticstudydata
    export get_gmap, get_alleles, get_chromosome, get_crossinfo, get_crosstype 
    export get_geno, get_genotype, get_genotranspose, get_pmap
    export get_phenocovar, get_pheno, get_isxchar
    
    #########
    # Plots #
    #########
    include("plots/plots_utils.jl")
    export gmap2df, pmap2df

    include("plots/plots_qtl.jl")
    export plot_QTL 

    include("plots/plots_manhattan.jl")
    export plot_manhattan

    include("plots/plots_eqtl.jl")
    export plot_eQTL 

end
