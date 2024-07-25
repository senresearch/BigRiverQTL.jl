# BigRiverQTL.jl



*A Statistical Toolbox for QTL Analysis*

`BigRiverQTL.jl` is a user-friendly Julia package that supports efficient and interpretable quantitative trait locus (QTL) analysis. This comprehensive toolbox encompasses three core components tailored to streamline the entire QTL analysis workflow: preprocessing, genomic scanning, and result visualization.

- **Preprocessing:** The preprocessing functions are designed to seamlessly import and convert genomic data into an efficient and memory-conservative format. This component also offers function capabilities for quickly calculating kinship matrices, ensuring data readiness for subsequent analysis phases.

- **Genomic Scanning:** `BigRiverQTL.jl` provides advanced genomic scanning capabilities through `BulkLMM.jl` for swift single-trait scans, which surpass other methods in terms of computational speed. For analyses involving multiple traits, the package employs `FlxQTL.jl`, a cutting-edge approach that detects complex trait interrelations.

- **Result Visualization:** The third component of `BigRiverQTL.jl` enriches the analytical experience by offering plotting tools designed to illustrate the outcomes of genomic scans. The plotting functions are useful for interpreting the results and aid in the derivation of meaningful conclusions


## Installation
To install `BigRiverQTL.jl`, you can use Julia's package manager. Here is the command:

```julia
using Pkg
Pkg.add("BigRiverQTL")
```



## Contribution
Contributions to BigRiverQTLPlots.jl are welcome and appreciated. If you'd like to contribute, please fork the repository and make changes as you'd like. If you have any questions or issues, feel free to open an 