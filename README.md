# Crisp
Crisp I: LCMS preprocessing, processing and multivariate statistical analysis

## Version history
2017-07-02 1st commit: R Scripts and Data files for preliminary analysis added. Readme.md added.

## Data files (in logical sequence order)
- Data5.RData (XCMS object with LC-MS data)
- injTable.RData (Injection table data)
- trtData.Rdata (Decoded feature matrix after within/between batch correction)
- ML_BigModels.RData (Multilevel analyses)
- permutation.RData (Permutation analyses H0 distributions)

## Script files (in logical sequence order)
- clustFuncBare2.r (Functions required for cluster-based within batch intensity drift correction)
- Drift1.R (Within batch corrections)
- Drift2.R (Between batch correction/normalisation)
- MLPLSVS_func.R (Functions required for multilevel analysis)
- ML-PLS-VS1.R (Multilevel multivariate statistical analysis incl permutations)
- ML-PLS-VS1.R (Data exctraction)
