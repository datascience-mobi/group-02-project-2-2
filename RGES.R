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
cds.treated <- newCountDataSet(countData=treated, condition=c(rep("treated",819))) 
cds.treated <- estimateSizeFactors(cds.treated)
cds.treated <- estimateDispersions(cds.treated, method="per-condition")

mode(untreated) <- "integer"
cds.untreated <- newCountDataSet(countData = untreated, condition = c(rep("untreated",819)))
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

p.values <- nbinomTestForMatrices(counts(cds.treated), counts(cds.untreated), sizeFactors(cds.treated), sizeFactors(cds.untreated), treated.dispersions[,1], untreated.dispersions[,1] )
names(p.values) <- row.names(counts(cds.treated))

# anderer Ansatz
mode(treated) <- "integer"
mode(untreated) <- "integer"
treated.untreated <- cbind(treated,untreated)
cds = newCountDataSet(countData = treated.untreated, conditions = c(rep("treated",819),rep("untreated",819)))
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
str( fitInfo(cds) )
plotDispEsts(cds)
fData(cds)
nbinom.treated.untreated = nbinomTest(cds, "treated", "untreated")
plotMA(nbinom.treated.untreated)
hist(nbinom.treated.untreated$pval, breaks=100, col="skyblue", border="slateblue", main="p-values nbinom")

