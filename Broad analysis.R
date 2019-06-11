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
colorvector = c("firebrick", "forestgreen", "blue", "orange", "black", "lightblue", "pink", "violet", "grey", "lightgreen", "navy", "gold", "yellow", "red", "peru")
drugnames <- drug
rownames(drugnames) = drugnames$Drug
drugcolor <- cbind(colorvector, rownames(drugnames))
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
boxplot(treated, xlab = "Samples", horizontal = F, border=drug.color.vector, main = "Gene Expression Treated")


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

#FC
log2FC.treated.untreated <- log2(treated.scaled/untreated.scaled)
is.nan.data.frame <- function(x)      #NaN durch 0 ersetzen
  do.call(cbind, lapply(x, is.nan))
log2FC.treated.untreated[is.nan(log2FC.treated.untreated)] <- 0
plot(density(log2FC.treated.untreated), main = "Log2 FC Treated/Untreated")

#PCA
PCA.FC <- prcomp(log2FC.treated.untreated, center=F , scale.=F)
plot(PCA.FC, type ="lines")
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], xlab = "PC1", ylab = "PC2", pch=19, main = "PCA Treated/Untreated")
# wie viel Varianz wird durch components erklaert?
Varianz.PCA=PCA.FC$sdev^2

# color: chemo/targeted
# 2 versions
# 1. sapply
#spalte ergänzen drug chemo
targeted.chemo <- c("targeted", "targeted", "chemo", "chemo", "chemo", "chemo", "chemo", "targeted", "targeted", "chemo", "chemo", "chemo", "chemo", "targeted", "chemo")
drug.added <- cbind(drug, targeted.chemo)
chemocolor <- ifelse(drug.added$targeted.chemo=="chemo", "firebrick", "forestgreen")
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
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=color.chemo, xlab = "PC1", ylab = "PC2", pch=19, main = "PCA Chemotherapy/Targeted")
plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=color.chemo, xlab = "PC3", ylab = "PC4", pch = 19, main = "PCA Chemotherapy/Targeted")

#plot colored PCA
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=drug.color.vector, xlab = "PC1", ylab = "PC2", pch=19, main ="PCA Drugs")
plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=drug.color.vector, xlab = "PC3", ylab = "PC4", pch=19, main = "PCA Drugs")

#color: tyrosine kinase inhibitor
# 2 versions
# 1. sapply
tyrosincolor <- ifelse(drug$Mechanism=="Tyrosine kinase inhibitor", "brown2", "darkolivegreen4")
tyrosincolordrugs <- cbind(tyrosincolor, rownames(drugnames))
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
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=color.tyrosin, xlab = "PC1", ylab = "PC2", pch=19, main = "PCA Tyrosin Kinase Inhibitor")
plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=color.tyrosin, xlab = "PC1", ylab = "PC2", pch=19, main = "PCA Tyrosin Kinase Inhibitor")





