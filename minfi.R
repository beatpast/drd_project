if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylation450kmanifest")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(minfi)
library(IlluminaHumanMethylation450kmanifest)