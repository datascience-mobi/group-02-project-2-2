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

#color drug: all drugs
# 2 versions
# 1. sapply
colorvector = c("firebrick", "forestgreen", "blue", "orange", "black", "lightblue", "pink", "violet", "grey", "lightgreen", "darkblue", "gold", "yellow", "red", "peru")
drugnames <- drug
rownames(drugnames) = drugnames$Drug
drugnames.ordered <- drugnames[order(drugnames$Drug),]
drugcolor <- cbind(colorvector, rownames(drugnames.ordered))
drugcolorvector <- sapply(rownames(meta), function(x){
  unname(drugcolor[meta[x, 3]], force = FALSE)
})
# 2. while loop
#linking the drugs to colors and creating a vector with the colors for the drugs to color the PCA
list.drug = list("5-Azacytidine", "bortezomib", "cisplatin","dasatinib","doxorubicin","erlotinib","geldanamycin","gemcitibine","lapatinib","paclitaxel","sirolimus","sorafenib","sunitinib","topotecan","vorinostat")

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

drug.color.vector = c(do.call("cbind",drug.color))

#general overview
#colored boxplot without scaling
boxplot(treated, xlab = "Samples", horizontal = F, border=drug.color.vector, main = "Gene Expression Treated", ylab= "Gene expression")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.6  , legend=rownames(drugnames.ordered), fill=colorvector, horiz=FALSE, cex=0.8, bty = "n")

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

#scale
basal.scaled <- scale(basalexp)
treated.scaled <- scale(treated)
untreated.scaled <- scale(untreated)

#colored boxplot with scaling
boxplot(treated.scaled, xlab="Samples", horizontal =F, border=drug.color.vector, main = "Gene Expression Treated Scaled")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.6  , legend=rownames(drugnames.ordered), fill=colorvector, horiz=FALSE, cex=0.8, bty = "n")

#FC
log2FC.treated.untreated <- log2(treated.scaled/untreated.scaled)
is.nan.data.frame <- function(x)      #NaN durch 0 ersetzen
  do.call(cbind, lapply(x, is.nan))
log2FC.treated.untreated[is.nan(log2FC.treated.untreated)] <- 0
plot(density(log2FC.treated.untreated), xlab= "log2 fold change values", main = "Density log2 fold change treated/untreated")

#PCA
PCA.FC <- prcomp(log2FC.treated.untreated, center=F , scale=F)
plot(PCA.FC, type ="lines", xlab ="principal component")
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], xlab = "PC1", ylab = "PC2", pch=19, main = "PCA for log2 fold change treated/untreated")
plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], xlab = "PC3", ylab = "PC4", pch=19, main = "PCA for log2 fold change treated/untreated")
# wie viel Varianz wird durch components erklaert?
Varianz.PCA=PCA.FC$sdev^2

# color: chemo/targeted
# 2 versions
# 1. sapply
#spalte ergänzen drug chemo
targeted.chemo <- c("targeted", "targeted", "chemo", "chemo", "chemo", "chemo", "chemo", "targeted", "targeted", "chemo", "chemo", "chemo", "chemo", "targeted", "chemo")
drug.added <- cbind(drug, targeted.chemo)
drug.added.ordered <- drug.added[order(drug.added$Drug),]
chemocolor <- ifelse(drug.added.ordered$targeted.chemo=="chemo", "firebrick", "forestgreen")
chemocolordrugs <- cbind(chemocolor, rownames(drugnames))
chemocolorvector <- sapply(rownames(meta), function(x){
       unname(chemocolordrugs[meta[x, 3]], force = FALSE)
   })
# 2. while loop
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
 

FC.named <- log2FC.treated.untreated
colnames(FC.named)[chemo] <- "chemo"
color.chemo <- ifelse(colnames(FC.named)=="chemo", "firebrick","forestgreen")

plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=color.chemo, xlab = "PC1", ylab = "PC2", pch=19, main = "PCA for log2 fold change treated/untreated")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.7  , legend=c("chemotherapy agent","targeted therapie"), fill=c("firebrick", "forestgreen"), horiz=FALSE, cex=0.8, bty = "n")

plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=color.chemo, xlab = "PC3", ylab = "PC4", pch = 19, main = "PCA for log2 fold change treated/untreated")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.7  , legend=c("chemotherapy agent","targeted therapie"), fill=c("firebrick", "forestgreen"), horiz=FALSE, cex=0.8, bty = "n")

#plot colored PCA
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=drug.color.vector, xlab = "PC1", ylab = "PC2", pch=19, main ="PCA for log2 fold change treated/untreated")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.6  , legend=rownames(drugnames.ordered), fill=colorvector, horiz=FALSE, cex=0.8, bty = "n")

plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=drug.color.vector, xlab = "PC3", ylab = "PC4", pch=19, main = "PCA for log2 fold change treated/untreated")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.6  , legend=rownames(drugnames.ordered), fill=colorvector, horiz=FALSE, cex=0.8, bty = "n")

#color: tyrosine kinase inhibitor
# 2 versions
# 1. sapply
tyrosincolor <- ifelse(drugnames.ordered$Mechanism=="Tyrosine kinase inhibitor", "brown2", "darkolivegreen4")
tyrosincolordrugs <- cbind(tyrosincolor, rownames(drugnames.ordered))
tyrosincolorvector <- sapply(rownames(meta), function(x){
          unname(tyrosincolor[meta[x, 3]], force = FALSE)
      })
# 2. while loop
list.tyrosin = list("dasatinib", "sunitinib", "lapatinib", "sorafenib")
tyrosin = c()
i=1
j=1
while(i<16)
{
  while(j<820)
  {
    if(isTRUE(meta[j,3]== list.chemo[i]))
    {tyrosin= c(tyrosin,j)
    }
    j = j +1
  }
  j= 1
  i=i+1
}

#plot tyrosin kinase PCA
colnames(FC.named)[tyrosin] <- "tyrosin"
color.tyrosin <- ifelse(colnames(FC.named)=="tyrosin", "brown2","darkolivegreen4")

plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=tyrosincolorvector, xlab = "PC1", ylab = "PC2", pch=19, main = "PCA Tyrosin Kinase Inhibitor")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.7  , legend=c("tyrosin kinase inhibitor","other"), fill=c("brown2", "darkolivegreen4"), horiz=FALSE, cex=0.8, bty = "n")

plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=tyrosincolorvector, xlab = "PC3", ylab = "PC4", pch=19, main = "PCA Tyrosin Kinase Inhibitor")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.7  , legend=c("tyrosin kinase inhibitor","other"), fill=c("brown2", "darkolivegreen4"), horiz=FALSE, cex=0.8, bty = "n")

