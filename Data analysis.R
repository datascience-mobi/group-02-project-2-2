library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)

meta = read.delim(paste0(wd, "/data/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t") 
untreated = readRDS(paste0(wd, "/data/NCI_TPW_gep_untreated.RDS"))
treated = readRDS(paste0(wd, "/data/NCI_TPW_gep_treated.RDS"))
ic50 = readRDS(paste0(wd, "/data/NegLogGI50.RDS"))
basalexp = readRDS(paste0(wd, "/data/CCLE_basalexpression.RDS"))
copynumber = readRDS(paste0(wd, "/data/CCLE_copynumber.RDS"))
mutations = readRDS(paste0(wd, "/data/CCLE_mutations.RDS"))
cellline = read.delim(paste0(wd, "/data/cellline_annotation.tsv"), header = TRUE, sep = "\t")
drug = read.delim(paste0(wd, "/data/drug_annotation.tsv"), header = TRUE, sep = "\t")

#general overview
#colored boxplot without scaling
boxplot(treated, xlab = "treated", horizontal = F, border=cb)
#colored boxplot with scaling
boxplot(treated.scaled, xlab="treated scaled", horizontal =F, border=cb)



#loop for NA values
list.na = list("treated"=treated, "untreated"=untreated, "mutations"=mutations, "basalexp" = basalexp, "cellline"=cellline, "copynumber"=copynumber, "ic50"=ic50, "meta"=meta)
myfunc = function(x){sum(is.na(x))}

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


# Spalten die chemo Medikamente enthalten
list.drug = list("5-Azacytidine", "bortezomib", "cisplatin","dasatinib","doxorubicin","erlotinib","geldanamycin","gemcitibine","lapatinib","paclitaxel","sirolimus","sorafenib","sunitinib","topotecan","vorinostat")
list.chemo = list("cisplatin","dasatinib","doxorubicin","geldanamycin","gemcitibine","lapatinib","paclitaxel","sorafenib","sunitinib","topotecan")
chemo = c()
i=1
j=1
while(i<16)
{
  while(j<820)
  {
       if(isTRUE(meta[j,3]== list.chemo[i]))
       {chemo= c(chemo,j)
         }
        j = j +1
  }
  j= 1
  i=i+1
}
 

# Medikamente einf�rben
FC.named <- log2FC.treated.untreated
colnames(FC.named)[chemo] <- "chemo"
color.chemo <- ifelse(colnames(FC.named)=="chemo", "firebrick","forestgreen")
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=color.chemo, xlab = "PC1", ylab = "PC2")

#color every drug PCA
FC.color <- log2FC.treated.untreated
#linking the drugs to colors and creating a vector with the colors for the drugs to color the PCA
list.drug = list("5-Azacytidine", "bortezomib", "cisplatin","dasatinib","doxorubicin","erlotinib","geldanamycin","gemcitabine","lapatinib","paclitaxel","sirolimus","sorafenib","sunitinib","topotecan","vorinostat")
list.colors = list("firebrick", "forestgreen", "blue", "orange", "black", "lightblue", "pink", "violet", "grey", "lightgreen", "darkblue", "gold", "yellow", "red", "peru", "gold")

i=1
j=1
a=1
drug.color = c()
while(i<16)
{
  while(j<820)
  {
    if(isTRUE(meta[j,3]== list.drug[i]))
    {drug.color = c(drug.color,list.colors[a])
    }
    j = j +1
  }
  j= 1
  i=i+1
  a=a+1
}

#creating color vector / convert the list drug.color into a vector
drug.color.vector = c(do.call("cbind",drug.color))
#plot colored PCA
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=cb, xlab = "PC1", ylab = "PC2", pch=19)
plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=cb, xlab = "PC3", ylab = "PC4", pch=19)

#color tyrosine kinase inhibitor
FC.TKI <- log2FC.treated.untreated
tyrosinKI <- c(dasatinib, sunitinib, lapatinib, sorafenib)
colnames(FC.TKI)[tyrosinKI]<- "Tyrosine Kinase Inhibitor"
color.TKI <- ifelse(colnames(FC.TKI)=="Tyrosine Kinase Inhibitor", "firebrick", "forestgreen")
#plot
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=drug.color.vector, xlab = "PC1", ylab = "PC2", pch=19)
plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=drug.color.vector, xlab = "PC3", ylab = "PC4", pch=19)


#Specific Analysisbo
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
upordown = t(upordown)
upordown.biomarker = upordown[ , c("SUPT3H","BBS9", "ANKRA2", "GNBL1", "KDM4C", "DDIT3", "TTC28","NBEA", "PCCA","PPP1R15A", "FTO", "PLK1", "VPS8", "POLR3B", "MAPKAP1", "STAG1", "TBCD", "C11orf49", "ANKS1A", "COMMD10")]
upordown.biomarker = t(upordown.biomarker)
heatmap(upordown.biomarker, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))

#biomarker Kriterium 2
biomarker.kriterium2 = log2FC.treated.untreated
biomarker.kriterium2[biomarker.kriterium2<0]=-1
biomarker.kriterium2[biomarker.kriterium2>0]=1
biomarker.kriterum2.sum = sapply(1:13299 , function(k) {sum(biomarker.kriterium2[k,]==-1)})
biomarker.kriterium2.sum.names =cbind(rownames(treated), biomarker.kriterum2.sum)
biomarker.kriterium2.sum.names.klein = biomarker.kriterium2.sum.names[,2]<205  #25%
biomarker.kriterium2.sum.names.gross = biomarker.kriterium2.sum.names[,2]>614
komplett.kriterium2 = cbind(biomarker.kriterium2.sum.names , biomarker.kriterium2.sum.names.klein , biomarker.kriterium2.sum.names.gross)
komplett.kriterium2 = komplett.kriterium2[-which(komplett.kriterium2[,4]=="FALSE" & komplett.kriterium2[,3]=="FALSE"),]
