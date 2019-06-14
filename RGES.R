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

### drug_signature Doris Versuch
mode(treated) <- "integer"
cds.treated <- newCountDataSet(countData=treated, condition=c(rep(1,819))) 
cds.treated <- estimateSizeFactors(cds.treated)
cds.treated <- estimateDispersions(cds.treated, method="per-condition")

mode(untreated) <- "integer"
cds.untreated <- newCountDataSet(countData = untreated, condition = c(rep(1,819)))
cds.untreated <- estimateSizeFactors(cds.untreated)
cds.untreated <- estimateDispersions(cds.untreated, method="per-condition")
# glaub brauchen wir eigentlich nicht von hier:
cols.treated <- conditions(cds.treated) == "A"
cols.untreated <- conditions(cds.untreated) == "B"

bmvA <- getBaseMeansAndVariances( counts(cds.treated)[,cols.treated], sizeFactors(cds.treated)[cols.treated] )
bmvB <- getBaseMeansAndVariances( counts(cds.untreated)[,cols.untreated], sizeFactors(cds.untreated)[cols.untreated] )
# bis hier.
str( fitInfo( cds.treated ) )
treated.dispersions <- fData(cds.treated)
str( fitInfo( cds.untreated ) )
untreated.dispersions <- fData(cds.untreated)

pvalues <- nbinomTestForMatrices(counts(cds.treated), counts(cds.untreated), sizeFactors(cds.treated), sizeFactors(cds.untreated), treated.dispersions[,1], untreated.dispersions[,1] )
names(pvalues) <- row.names(counts(cds.treated))
