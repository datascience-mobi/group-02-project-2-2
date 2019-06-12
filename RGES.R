### package for nbinomTest
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq")

###data
library(rstudioapi)
wd = dirname(rstudioapi::getSourceEditorContext()$path)
meta = read.delim(paste0(wd, "/data/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t") 
untreated = readRDS(paste0(wd, "/data/NCI_TPW_gep_untreated.RDS"))
treated = readRDS(paste0(wd, "/data/NCI_TPW_gep_treated.RDS"))
basalexp = readRDS(paste0(wd, "/data/CCLE_basalexpression.RDS"))
library(DESeq)

### drug_signature
treated=estimateSizeFactorsForMatrix(treated)
treated=estimateDispersions(treated, method= "per-condition") #problem
untreated=estimateSizeFactorsForMatrix(untreated)
untreated=estimateDispersions(untreated, method="per-condition") #problem
nbinomTestForMatrices(treated, untreated, ...)
