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



###basal anpassen - jojos versuch - die neue heisst jetzt basal.fitted.untreated
# nurnoch die Gene dalassen welche in basal und untreated sind
basal.fit = subset(basalexp, rownames(basalexp) %in% rownames(untreated))
untreated.fit = subset(untreated, rownames(untreated) %in% rownames(basalexp))

meta.matrix <- as.matrix(meta)
dataset.basal <- as.matrix(basal.fit)

new.basal.names <- as.character(meta.matrix[1:819,2])
output.dataset <- sapply(seq_along(new.basal.names), function(a) {
  name_picker <- new.basal.names[a]
  out <- dataset.basal[,which(colnames(dataset.basal) == name_picker)]
  return(out)
})

#als Matrix umformatieren und umnennen damit wir es erkennen
basal.fitted.untreated <- matrix(unlist(output.dataset), nrow = 11461, ncol = 819, byrow=FALSE, dimnames = NULL)
colnames(basal.fitted.untreated) <- new.basal.names

# rownames: gene einfügen
rownames(basal.fitted.untreated)=rownames(basal.fit)

#signature genes
#die Kriterien von Bin Chen schmeißen bei uns alle Gene raus
#for the RGESexample code they use 978 genes so we aim for the same number of genes
#log2 Kriterium ganz raus weil dafür sind unsere Werte viel zu klein
dz_signature <- subset(nbinom.treated.untreated, !is.na(padj) & !is.na(id) & id !='?' & padj < 0.5  & abs(log2FoldChange) != Inf )
dim(dz_signature)
gene.list = c(dz_signature[,1])
  
#log2FC 
#natürlich auch noch die für untreated/basal das muss dann signature heißen
drug_signature <- log2(treated/untreated)
is.nan.data.frame <- function(x)     
  do.call(cbind, lapply(x, is.nan))
drug_signature[is.nan(drug_signature)] <- 0

#only keep FC values for dz_signature genes 
drug_signature = subset(drug_signature, rownames(drug_signature) %in% gene.list)
signature = subset(signature, rownames(signature) %in% gene.list)
