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

#general overview
#colored boxplot without scaling
boxplot(treated, xlab = "treated", horizontal = F, border=cb)
#colored boxplot with scaling
boxplot(treated.scaled, xlab="treated scaled", horizontal =F, border=cb)

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

#scaled, FC, PCA
is.nan.data.frame <- function(x)      #NaN durch 0 ersetzen
  do.call(cbind, lapply(x, is.nan))
log2FC.treated.untreated[is.nan(log2FC.treated.untreated)] <- 0
plot(density(log2FC.treated.untreated))
PCA.FC <- prcomp(log2FC.treated.untreated, center=F , scale.=F)
plot(PCA.FC, type ="lines")
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], xlab = "PC1", ylab = "PC2")

# wie viel Varianz wird durch components erklärt?
PCA.FC$sdev^2

# Medikamente einfï¿½rben
colnames(FC.named)[chemo] <- "chemo"
color.chemo <- ifelse(colnames(FC.named)=="chemo", "firebrick","forestgreen")
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=color.chemo, xlab = "PC1", ylab = "PC2")
#color every drug PCA
FC.color <- log2FC.treated.untreated
#creating vectors for every drug and rename the log2FC matrix
azacytidine <- c(1:39)
colnames(FC.color)[azacytidine] <- "azacytidine"
bortezomib <-c(40:94)
colnames(FC.color)[bortezomib] <- "bortezomib"
cisplatin <- c(95:150)
colnames(FC.color)[cisplatin] <- "cisplatin"
dasatinib <- c(150:198)
colnames(FC.color)[dasatinib] <- "dasatinib"
doxorubicin <- c(199:248)
colnames(FC.color)[doxorubicin] <- "doxorubicin"
erlotinib <- c(249:307)
colnames(FC.color)[erlotinib] <- "erlotinib"
geldanamycin <- c(308:364)
colnames(FC.color)[geldanamycin] <- "geldanamycin"
gemcitabine <- c(365:420)
colnames(FC.color)[gemcitabine] <- "gemcitabine"
lapatinib <- c(421:474)
colnames(FC.color)[lapatinib] <- "lapatinib"
paclitaxel <- c(475:533)
colnames(FC.color)[paclitaxel] <- "paclitaxel"
sirolimus <- c(534:589)
colnames(FC.color)[sirolimus] <- "sirolimus"
sorafenib <- c(590:646)
colnames(FC.color)[sorafenib] <- "sorafenib"
sunitinib <- c(647:702)
colnames(FC.color)[sunitinib] <- "sunitinib"
topotecan <- c(703:760)
colnames(FC.color)[topotecan] <- "topotecan"
vorinostat <- c(761:819)
colnames(FC.color)[vorinostat] <- "vorinostat"
#creating color vector
color.aza <- ifelse(colnames(FC.color[,c(1:94)])=="azacytidine", "firebrick", "forestgreen")
color.cis <- ifelse(colnames(FC.color[,c(95:198)])=="cisplatin", "blue","orange")
color.dox <- ifelse(colnames(FC.color[,c(199:307)])=="doxorubicin", "black", "lightblue")
color.gel <- ifelse(colnames(FC.color[,c(308:420)])=="geldanamycin", "pink", "violet")
color.lap <- ifelse(colnames(FC.color[,c(421:533)])=="lapatinib", "grey", "lightgreen")
color.sir<- ifelse(colnames(FC.color[,c(534:646)])=="sirolimus", "darkblue", "gold")
color.sun <- ifelse(colnames(FC.color[,c(647:760)])=="sunitinib", "yellow", "red")
color.vor <- ifelse(colnames(FC.color[,c(761:819)])=="vorinostat", "peru", "silver")
cb <- c(color.aza, color.cis, color.dox, color.gel, color.lap, color.sir, color.sun, color.vor)
#plot colored PCA
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=cb, xlab = "PC1", ylab = "PC2", pch=19)
plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=cb, xlab = "PC3", ylab = "PC4", pch=19)

#color tyrosine kinase inhibitor
FC.TKI <- log2FC.treated.untreated
tyrosinKI <- c(dasatinib, sunitinib, lapatinib, sorafenib)
colnames(FC.TKI)[tyrosinKI]<- "Tyrosine Kinase Inhibitor"
color.TKI <- ifelse(colnames(FC.TKI)=="Tyrosine Kinase Inhibitor", "firebrick", "forestgreen")
#plot
plot(PCA.FC$rotation[,1], PCA.FC$rotation[,2], col=color.TKI, xlab="PC1", ylab="PC2", pch=19)
plot(PCA.FC$rotation[,3], PCA.FC$rotation[,4], col=color.TKI, xlab="PC3", ylab="PC4", pch=19)
