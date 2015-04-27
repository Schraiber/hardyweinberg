# hardweinberg

This set of R functions will perform maximum likelihood estimation of inbreeding coefficients in small sample sizes by fitting the exact distribution of the number of heterozygotes. Further, it can perform composite likelihood estimation of inbreeding coefficients in a single individual.

# Using the scripts

## Installation

Calling `source("HWtest.r")` to load the scripts. The only required package is `EMT`.

## Data input format

Data is required to be in a format similar, but not identical to, a VCF file. The required columns are: 1) chromosome, 2) position, 3) number of heterozygotes, 4) allele frequency, 5-end) a 1 for each individual who is a heterozygote, and a 0 for each individual who is homozygous (either reference or alternative). For example, the following data corresponds to the 9 YRI individuals in the Complete Genomics dataset:

```
 V1     V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13
  1 564654  0  8  0  0  0  0  0   0   0   0   0
  1 564936  0  2  0  0  0  0  0   0   0   0   0
  1 565317  0  2  0  0  0  0  0   0   0   0   0
  1 565371  1  1  0  0  0  0  0   1   0   0   0
  1 713698  2  2  0  1  0  0  0   0   1   0   0
  1 721119  1  1  0  0  0  1  0   0   0   0   0

```

## Batch inbreeding coefficient

To compute the liklelihood of the data assuming  a single inbreeding coefficient for all individuals in the dataset, use the function `LL\_combined`.
