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



#gereral overview
##creating a color vektor

# 1. sapply
colorvector = c("firebrick", "forestgreen", "blue", "orange", "black", "lightblue", "pink", "violet", "grey", "lightgreen", "darkblue", "gold", "yellow", "red", "peru")
drugnames <- drug
rownames(drugnames) = drugnames$Drug
drugnames.ordered <- drugnames[order(drugnames$Drug),]
drugcolor <- cbind(colorvector, rownames(drugnames.ordered))
drugcolorvector <- sapply(rownames(meta), function(x){
  unname(drugcolor[meta[x, 3]], force = FALSE)
})


##boxplot with the samples from treated, color every drug in the boxplot without scaling
boxplot(treated, xlab = "Samples", horizontal = F, border=drugcolorvector, main = "Gene Expression Treated", ylab= "Gene expression")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.3  , legend=rownames(drugnames.ordered), fill=colorvector, horiz=FALSE, cex=0.8, bty = "n")


#scale
basal.scaled <- scale(basalexp)
treated.scaled <- scale(treated)
untreated.scaled <- scale(untreated)

##colored boxplot with scaling
boxplot(treated.scaled, xlab="Samples", horizontal =F, border=drugcolorvector, main = "Gene Expression Treated Scaled")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.2  , legend=rownames(drugnames.ordered), fill=colorvector, horiz=FALSE, cex=0.8, bty = "n")

##compute the FC
log2FC.treated.untreated <- log2(treated.scaled/untreated.scaled)
is.nan.data.frame <- function(x)      #NaN durch 0 ersetzen
  do.call(cbind, lapply(x, is.nan))
log2FC.treated.untreated[is.nan(log2FC.treated.untreated)] <- 0
plot(density(log2FC.treated.untreated), xlab= "log2 fold change values", main = "Density log2 fold change treated/untreated")

#PCA
PCA.FC <- prcomp(log2FC.treated.untreated, center=F , scale=F)
plot(PCA.FC, type ="lines", main = "Elbow plot of PCA")
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], xlab = "PC1", ylab = "PC2", pch=19, main = "PCA for log2 fold change treated/untreated")
plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], xlab = "PC3", ylab = "PC4", pch=19, main = "PCA for log2 fold change treated/untreated")
# wie viel Varianz wird durch components erklaert?
Varianz.PCA=PCA.FC$sdev^2

# color the pca: chemo/targeted
# 2 versions

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


##plotting pca - color chemo targeted

plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=color.chemo, xlab = "PC1", ylab = "PC2", pch=19, main = "PCA for log2 fold change treated/untreated")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.33  , legend=c("chemotherapy","targeted therapy"), fill=c("firebrick", "forestgreen"), horiz=FALSE, cex=0.8, bty = "n")

plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=color.chemo, xlab = "PC3", ylab = "PC4", pch = 19, main = "PCA for log2 fold change treated/untreated")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.33  , legend=c("chemotherapy","targeted therapy"), fill=c("firebrick", "forestgreen"), horiz=FALSE, cex=0.8, bty = "n")


#plot colored PCA - color every drug

plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=drugcolorvector, xlab = "PC1", ylab = "PC2", pch=19, main ="PCA for log2 fold change treated/untreated")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.3  , legend=rownames(drugnames.ordered), fill=colorvector, horiz=FALSE, cex=0.8, bty = "n")

plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=drugcolorvector, xlab = "PC3", ylab = "PC4", pch=19, main = "PCA for log2 fold change treated/untreated")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.3  , legend=rownames(drugnames.ordered), fill=colorvector, horiz=FALSE, cex=0.8, bty = "n")

#plot colored PCA - color tyrosine kinase inhibitor
# 2 versions
# 1. sapply
tyrosincolor <- ifelse(drugnames.ordered$Mechanism=="Tyrosine kinase inhibitor", "brown2", "darkolivegreen4")
tyrosincolordrugs <- cbind(tyrosincolor, rownames(drugnames.ordered))
tyrosincolorvector <- sapply(rownames(meta), function(x){
          unname(tyrosincolor[meta[x, 3]], force = FALSE)
      })


#plot

plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], col=tyrosincolorvector, xlab = "PC1", ylab = "PC2", pch=19, main = "PCA Tyrosin Kinase Inhibitor")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.27  , legend=c("TK inhibitor","other"), fill=c("brown2", "darkolivegreen4"), horiz=FALSE, cex=0.8, bty = "n")

plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], col=tyrosincolorvector, xlab = "PC3", ylab = "PC4", pch=19, main = "PCA Tyrosin Kinase Inhibitor")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend( x= "right" ,inset=-0.27  , legend=c("TK inhibitor","other"), fill=c("brown2", "darkolivegreen4"), horiz=FALSE, cex=0.8, bty = "n")





######PCA nach JUli anfärben
library(ggplot2)
library(viridis)

pca = prcomp(log2FC.treated.untreated)
summary(pca)

meta_neu = meta[1:819,]
meta_neu = as.data.frame(meta_neu)

#plot(pca$x[,2], pca$x[,3], col = meta_neu$drug)
plot(pca$rotation[,1], pca$rotation[,3], by=meta_neu$drug, pch = 19, cex = 0.8, cols = viridis(15) )




##################### BABY ##########
library(ggplot2)
library(reshape2)
##
melt = melt(treated)
drugcolor = vector()
x=1
while(x < 820){
  drugcolor = c(drugcolor, rep( as.character(meta[x,3]) , 13299))
  x=x+1
}
melt = cbind(melt ,drugcolor)

cols15 <- c("#d7ff2e","#2f0085","#02ee81","#f357ff",
            "#86ae00","#ff2f86","#008118","#ff443e",
            "#008b9d","#ffcd77","#1a001e","#f8ffb4",
            "#b6a7ff","#572a00","#c8fff7")

ggplot(melt, aes(x = Var2, y = value)) + labs(y='Gene') + theme_bw() +
  geom_boxplot(aes(fill = melt$drugcolor, color=melt$drugcolor) ,outlier.size = 0.1) +
  ggtitle("Gene Expression Treated") +
  scale_fill_manual(name= melt$drugcolor, values = cols15) + 
  scale_color_manual(name = melt$drugcolor , values = cols15)

