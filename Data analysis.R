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





#loop for NA values
list.na = list("treated"=treated, "untreated"=untreated, "mutations"=mutations, "basalexp" = basalexp, "cellline"=cellline, "copynumber"=copynumber, "ic50"=ic50, "meta"=meta)
myfunc = function(x){sum(is.na)}

i = 0
while(i<9)
{
  p=sapply(list.na[i],myfunc)
  print(p)
  i = i +1
}

mutations.removed = mutations[, -(12:13)]

#scaled, FC, PCA
basal.scaled <- scale(basalexp)
treated.scaled <- scale(treated)
untreated.scaled <- scale(untreated)
log2FC.treated.untreated <- log2(treated.scaled/untreated.scaled)

is.nan.data.frame <- function(x)      #NaN durch 0 ersetzen
  do.call(cbind, lapply(x, is.nan))
log2FC.treated.untreated[is.nan(log2FC.treated.untreated)] <- 0
plot(density(log2FC.treated.untreated))
PCA.FC <- prcomp(log2FC.treated.untreated, center=F , scale.=F)
plot(PCA.FC, type ="lines")
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], xlab = "PC1", ylab = "PC2")

# wie viel Varianz wird durch components erklaert?
Varianz.PCA=PCA.FC$sdev^2


# Medikamente einfaerben
chemo1 <- c(94:248)   #ueberpruefen
chemo2 <- c(307:533)  #ueberpruefen
chemo3 <- c(589:760)  #ueberpruefen
chemo <- c(chemo1, chemo2, chemo3)
FC.named <- log2FC.treated.untreated


# Medikamente einf�rben
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


#Specific Analysis
#Find Biomarker for cisplatin through FC
FC.cisplatin = log2FC.treated.untreated[,c(95:149)]
FC.cisplatin = t(FC.cisplatin)
mean.FC <- colMeans(FC.cisplatin, na.rm = TRUE, dims =1)
sorted.mean = sort(mean.FC, decreasing = TRUE)

#finding/visualizing most extreme FC values for cisplatin
lowest.FC <- sorted.mean[13290:13299]
par(mar = c(5, 5, 5, 5))
barplot(lowest.FC, horiz = TRUE, xlim = c(-1.3,0), main= "lowest log2 FC-values for cisplatin", xlab= "mean log2FC values in different celllines", col= "firebrick", names.arg=c("FTO", "PLK1", "VPS8", "POLR38", "MAPKAP1", "STAG1", "TBCD", "C11orf49","ANKS1A", "COMMD10"), las=1, cex.names =0.8)

highest.FC <- sorted.mean[1:10]
par(mar = c(5, 5, 5, 5))
barplot(highest.FC, horiz = TRUE, xlim = c(0,1.8), main= "highest log2 FC-values for cisplatin", xlab= "mean log2FC values in different celllines", col= "lightgreen",las=1, cex.names =0.8)

#upordownregulation
upordown = (treated.scaled/untreated.scaled)
upordown = upordown[,c(95:149)]
upordown= t(upordown)
upordown.biomarker = upordown[,c(11372, 1053, 503, 1347, 4403, 5711, 2819, 12310, 7537, 7538, 7539, 8326, 9002, 4060, 8765, 12703, 8892, 6729, 11254, 6344, 11536, 1261, 525, 2332)]
heatmap(upordown.biomarker)

