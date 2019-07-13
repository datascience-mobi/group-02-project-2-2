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
##Boxplots
###boxplot with the samples from treated, color every drug in the boxplot without scaling

library(ggplot2)
library(reshape2)

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


##scale
basal.scaled <- scale(basalexp)
treated.scaled <- scale(treated)
untreated.scaled <- scale(untreated)

###colored boxplot with scaling

library(ggplot2)
library(reshape2)
##
melt.scale = melt(treated.scaled)
drugcolor = vector()
x=1
while(x < 820){
  drugcolor = c(drugcolor, rep( as.character(meta[x,3]) , 13299))
  x=x+1
}
melt.scale = cbind(melt.scale ,drugcolor)

cols15 <- c("#d7ff2e","#2f0085","#02ee81","#f357ff",
            "#86ae00","#ff2f86","#008118","#ff443e",
            "#008b9d","#ffcd77","#1a001e","#f8ffb4",
            "#b6a7ff","#572a00","#c8fff7")

ggplot(melt.scale, aes(x = Var2, y = value)) + labs(y='Gene') + theme_bw() +
  geom_boxplot(aes(fill = melt.scale$drugcolor, color=melt.scale$drugcolor) ,outlier.size = 0.1) +
  ggtitle("Gene Expression Treated") +
  scale_fill_manual(name= melt.scale$drugcolor, values = cols15) + 
  scale_color_manual(name = melt.scale$drugcolor , values = cols15)


##compute the FC values
log2FC.treated.untreated <- treated.scaled - untreated.scaled
is.nan.data.frame <- function(x)      #NaN durch 0 ersetzen
  do.call(cbind, lapply(x, is.nan))
log2FC.treated.untreated[is.nan(log2FC.treated.untreated)] <- 0
plot(density(log2FC.treated.untreated), xlab= "log2 fold change values", main = "Density log2 fold change treated/untreated")

##PCA

PCA.FC <- prcomp(log2FC.treated.untreated, center=F , scale=F)
plot(PCA.FC, type ="lines", main = "Elbow plot of PCA")
plot(PCA.FC$rotation[, 1], PCA.FC$rotation[, 2], xlab = "PC1", ylab = "PC2", pch=19, main = "PCA for log2 fold change treated/untreated")
plot(PCA.FC$rotation[, 3], PCA.FC$rotation[, 4], xlab = "PC3", ylab = "PC4", pch=19, main = "PCA for log2 fold change treated/untreated")

### wie viel Varianz wird durch components erklaert?
Varianz.PCA=PCA.FC$sdev^2

###plotting PCA - color every drug

library(ggplot2)
library(viridis)

pca = prcomp(log2FC.treated.untreated)
summary(pca)

meta_neu = meta[1:819,]
meta_neu = as.data.frame(meta_neu)

pca_plot_drugs12 <- ggplot(as.data.frame(pca$rotation), aes(x= pca$rotation[,1], y = pca$rotation[,2])) +
  theme_bw(base_size = 7) +
  geom_point(aes(colour = factor(meta_neu$drug))) +
  scale_colour_viridis(option ="viridis", discrete = TRUE) +
  #mit BREWER: scale_fill_brewer(palette = "Dark2")
  ggtitle("Principal Component Analysis 1 & 2") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2")

pca_plot_drugs12

pca_plot_drugs34 <- ggplot(as.data.frame(pca$rotation), aes(x= pca$rotation[,3], y = pca$rotation[,4])) +
  theme_bw(base_size = 7) +
  geom_point(aes(colour = factor(meta_neu$drug))) +
  scale_colour_viridis(option ="viridis", discrete = TRUE) +
  #mit BREWER: scale_fill_brewer(palette = "Dark2")
  ggtitle("Principal Component Analysis 3 & 4") +
  xlab("Principal Component 3") +
  ylab("Principal Component 4")

pca_plot_drugs34

###plotting PCA - color chemo targeted
####an meta_neu Spalte anfügen, ob chemo oder targeted; undzwar, wenn in meta_neu spalte drug = drug.addad.ordered; dann füge
#### in meta in neue spalte das ein, was in drug.adde.ordered in targeted.chemo steht

targeted.chemo <- c("targeted", "targeted", "chemo", "chemo", "chemo", "chemo", "chemo", "targeted", "targeted", "chemo", "chemo", "chemo", "chemo", "targeted", "chemo")
drug.added <- cbind(drug, targeted.chemo)
drug.added.ordered <- drug.added[order(drug.added$Drug),]

chem.targ = c(rep(as.numeric(0),819))
meta_neu = cbind(meta_neu, chem.targ)

i=1
j=1

while(j<16)
{while(i<820)
{
  if(isTRUE(meta_neu[i,3]== drug.added.ordered[j,1]))
  {meta_neu[i,7] = as.character(drug.added.ordered[j,9])
  }
  i = i +1
}
  i= 1
  j=j+1}

###pca

pca_plot_chemtarg12 <- ggplot(as.data.frame(pca$rotation), aes(x= pca$rotation[,1], y = pca$rotation[,2])) +
  theme_bw(base_size = 7) +
  geom_point(aes(colour = factor(meta_neu$chem.targ))) +
  scale_colour_viridis(option ="viridis", discrete = TRUE) +
  #mit BREWER: scale_fill_brewer(palette = "Dark2")
  ggtitle("Principal Component Analysis 1&2") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2")

pca_plot_chemtarg12

pca_plot_chemtarg34 <- ggplot(as.data.frame(pca$rotation), aes(x= pca$rotation[,3], y = pca$rotation[,4])) +
  theme_bw(base_size = 7) +
  geom_point(aes(colour = factor(meta_neu$chem.targ))) +
  scale_colour_viridis(option ="viridis", discrete = TRUE) +
  #mit BREWER: scale_fill_brewer(palette = "Dark2")
  ggtitle("Principal Component Analysis 3&4") +
  xlab("Principal Component 3") +
  ylab("Principal Component 4")

pca_plot_chemtarg34

###plotting PCA - color tyrosine kinase inhibitor
###an meta_neu Spalte anfügen, ob TKI oder nicht; undzwar, wenn in meta_neu spalte drug = drug.added.ordered UND in
### drug.added.ordered Spalte Mechanism=TKI ; dann füge in meta in neue Spalte das ein was bei Mechanism steht

TKI = c(rep(as.numeric(0),819))
meta_neu = cbind(meta_neu, TKI)

i=1
j=1

while(j<16)
{while(i<820)
{
  if(isTRUE(meta_neu[i,3]== drug.added.ordered[j,1])
     & (drug.added.ordered[j,3] == "Tyrosine kinase inhibitor"))
  {meta_neu[i,8] = as.character(drug.added.ordered[j,3])
  }
  i = i +1
}
  i= 1
  j=j+1}

##pca

pca_plot_TKI12 <- ggplot(as.data.frame(pca$rotation), aes(x= pca$rotation[,1], y = pca$rotation[,2])) +
  theme_bw(base_size = 7) +
  geom_point(aes(colour = factor(meta_neu$TKI))) +
  scale_colour_viridis(option ="viridis", discrete = TRUE) +
  ggtitle("Principal Component Analysis 1&2") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2")

pca_plot_TKI12

pca_plot_TKI34 <- ggplot(as.data.frame(pca$rotation), aes(x= pca$rotation[,3], y = pca$rotation[,4])) +
  theme_bw(base_size = 7) +
  geom_point(aes(colour = factor(meta_neu$TKI))) +
  scale_colour_viridis(option ="viridis", discrete = TRUE) +
  ggtitle("Principal Component Analysis 3&4") +
  xlab("Principal Component 3") +
  ylab("Principal Component 4")

pca_plot_TKI34

