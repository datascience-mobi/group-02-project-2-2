library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)

meta = read.delim(paste0(wd, "/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t") 
untreated = readRDS(paste0(wd, "/NCI_TPW_gep_untreated.RDS"))
treated = readRDS(paste0(wd, "/NCI_TPW_gep_treated.RDS"))
ic50 = readRDS(paste0(wd, "/NegLogGI50.RDS"))
basalexp = readRDS(paste0(wd, "/CCLE_basalexpression.RDS"))
copynumber = readRDS(paste0(wd, "/CCLE_copynumber.RDS"))
mutations = readRDS(paste0(wd, "/CCLE_mutations.RDS"))
cellline = read.delim(paste0(wd, "/cellline_annotation.tsv"), header = TRUE, sep = "\t")
drug = read.delim(paste0(wd, "/drug_annotation.tsv"), header = TRUE, sep = "\t")

na.treated = apply(treated, 2, function(x) {sum(is.na(x))})
which(na.treated > 0)
na.untreated = apply(untreated, 2, function(x) {sum(is.na(x))})
which(na.untreated > 0)
na.mutations = apply(mutations, 2, function(x) {sum(is.na(x))})
which(na.mutations > 0)
na.basal=apply(basalexp, 2, function(x) {sum(is.na(x))})
which(na.basal > 0)
na.cellline=apply(cellline, 2, function(x) {sum(is.na(x))})
which(na.cellline>0)
na.copynumber=apply(copynumber, 2, function(x) {sum(is.na(x))})
which(na.copynumber > 0)
na.ic50=apply(ic50, 2, function(x) {sum(is.na(x))})
which(ic50 > 0)
na.meta=apply(meta, 2, function(x) {sum(is.na(x))})
which(na.meta >0)

mutations.removed = mutations[, -(12:13)]
basal.scaled <-- scale(basalexp)
treated.scaled <-- scale(treated)
untreated.scaled <-- scale(untreated)
log2FC.treated.untreated <-- log2(treated.scaled/untreated.scaled)
is.nan.data.frame <- function(x)      #NaN durch 0 ersetzen
  do.call(cbind, lapply(x, is.nan))
log2FC.treated.untreated[is.nan(log2FC.treated.untreated)] <- 0
plot(density(log2FC.treated.untreated))
PCA.FC <- prcomp(log2FC.treated.untreated, center=F , scale.=F)
plot(PCA.FC, type ="lines")
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], xlab = "PC1", ylab = "PC2")
# Medikamente einf�rben
