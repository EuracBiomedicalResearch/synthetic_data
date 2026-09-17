# Extensions to the HAPNEST workflow

This fork extends the original HAPNEST workflow with new functionality
for the generation and evaluation of large-scale synthetic
genotype–phenotype datasets, together with stability, compatibility,
and reproducibility improvements.

## New implementations

### GWAS Catalog integration
- Added automated integration of GWAS Catalog associations for phenotype
  simulation.
- This functionality supports the definition of phenotype-specific causal
  variants using associations derived from the GWAS Catalog.

### VCF export
- Added support for exporting generated synthetic genotype data in Variant
  Call Format (VCF).
- This enables synthetic genotype datasets to be distributed and analyzed
  using a standard format supported by commonly used genomics tools.

### Membership Inference Attack evaluation
- Added Membership Inference Attack (MIA) analyses for privacy-oriented
  evaluation of the generated synthetic genomic datasets.

## Stability, compatibility, and bug fixes

### Reproducible Julia environment
- Added `Manifest.toml` to preserve the exact Julia package versions used
  by the original HAPNEST environment.
- This prevents dependency resolution from installing incompatible package
  versions when rebuilding the environment.

### Compatibility improvements
- Resolved software, dependency, and execution-environment compatibility
  issues identified during workflow deployment and testing.
- Added support for both AVX2 and generic x86_64 PLINK2 binaries, with
  automatic runtime selection according to the CPU capabilities of the
  execution node. This enables execution on heterogeneous compute nodes
  with and without AVX2 support.

### PCA evaluation
- Replaced the KING-based PCA projection step with a PLINK2-based projection
  using `plink2 --score`. The original KING-based projection was found to
  be unreliable for large synthetic cohorts, causing segmentation faults
  with large sample counts, consistent with an internal 32-bit integer
  overflow.
- Added dynamic detection of the relevant columns in the PLINK2 PCA output
  to ensure compatibility across PLINK2 versions and configurations.
- Added normalization of projected PC scores to make them comparable with
  the reference PCA coordinates.
  
### GWAS numerical stability
- Added covariate variance standardization to `plink2 --glm` to prevent
  numerical stability failures observed for large synthetic cohorts.

### GWAS P-value handling
- Improved handling of extremely small GWAS P-values falling outside the
  standard Float64 range, preventing failures during downstream GWAS
  visualization.

## Scope

These modifications extend the functionality of the original HAPNEST
workflow while improving its stability, compatibility, and reproducibility
for the generation and evaluation of population-specific synthetic
genotype–phenotype datasets at large sample sizes.

For the original workflow description and general usage instructions,
see the main repository documentation.