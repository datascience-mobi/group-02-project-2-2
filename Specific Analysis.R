## packages 
library(rstudioapi)
library(gplots)
library(pheatmap)

## necessary data from the general analysis
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

#scale
basal.scaled <- scale(basalexp)
treated.scaled <- scale(treated)
untreated.scaled <- scale(untreated)

#FC
log2FC.treated.untreated <- log2(treated.scaled/untreated.scaled)
is.nan.data.frame <- function(x)      #NaN durch 0 ersetzen
  do.call(cbind, lapply(x, is.nan))
log2FC.treated.untreated[is.nan(log2FC.treated.untreated)] <- 0


##Specific analysis

## Step 1 (a)

#Find Biomarker for cisplatin through FC (criterium 1)
cisplatin.col=c()
j=1

while(j<820)
{
  if(isTRUE(meta[j,3]== "cisplatin"))
  {cisplatin.col= c(cisplatin.col,j)
  }
  j = j +1
}

FC.cisplatin = log2FC.treated.untreated[,cisplatin.col]
mean.FC.rows <- apply(FC.cisplatin, 1, mean)
sorted.mean = sort(mean.FC.rows, decreasing = TRUE)

#finding/visualizing most extreme FC values for cisplatin
lowest.FC <- sorted.mean[13280:13299]
par(mar = c(5, 5, 5, 5))
barplot(lowest.FC,
        horiz = TRUE,
        xlim = c(-1.3,0),
        main= "lowest log2 FC-values for cisplatin",
        xlab= "mean log2FC values in different celllines",
        col= "firebrick",
        las=1,
        cex.names =0.8)

highest.FC <- sorted.mean[1:20]
par(mar = c(5, 5, 5, 5))
barplot(highest.FC,
        horiz = TRUE,
        xlim = c(0,1.8),
        main= "highest log2 FC-values for cisplatin",
        xlab= "mean log2FC values in different celllines",
        col= "lightgreen",
        las=1,
        cex.names =0.8)


#matrix containing the biomarker found through the FC 
highest.FC = as.matrix(highest.FC) 
lowest.FC = as.matrix(lowest.FC)
biomarker1.FC = as.matrix(c(highest.FC, lowest.FC))
highest.names <- row.names(highest.FC)
lowest.names <- row.names(lowest.FC)
row.names(biomarker1.FC) <- c(highest.names, lowest.names)
biomarker1 <- c(highest.names, lowest.names)

#criterium 2
is.neg = FC.cisplatin<0
i =1
j=1
a=1
biomarker2.up = c()
biomarker2.down = c()
while(j<13300){
  while(i<56){
    if(is.neg[j,i]==TRUE)
     {a=a+1}
    i=i+1}
  if(a>49)
     {biomarker2.down= c(biomarker2.down, j)}
  if(a<6)
    {biomarker2.up= c(biomarker2.up, j)}
  a=1
  i=1
  j=j+1}

#das sind jetzt die nrow von Genen dich sich unter cisplatin
#in 50 von 55 Faellen in die gleiche Richtung veraendern
biomarker2 = row.names(FC.cisplatin[c(biomarker2.down, biomarker2.up),])

#do we find the same biomarkers for cisplatin with both criteria? 
i=1
j=1
a=1
double.biomarker = c()
while(i<41)
{
  while(j<375)
      {
        if(isTRUE(biomarker2[j] == biomarker1[i]))
           {double.biomarker = c(double.biomarker, biomarker1[a])
             }
            j = j +1
         }
          j= 1
          i=i+1
          a=a+1
}


## Step 1 (b)
# ttest to verify significance of the biomarker
#Creating seperated matrices containing cisplatin
treated.cisplatin <- treated.scaled[,grep ("cisplatin", colnames(treated.scaled))]
untreated.cisplatin <- untreated.scaled[,grep("cisplatin", colnames(treated.scaled))]

# 1. checking normality example BBS9
qqnorm(treated.cisplatin["BBS9", ], main = "BBS9")
qqline(treated.cisplatin["BBS9", ])

# 2. Welch two sample t-test
pvalues.welch <- sapply(double.biomarker, function(x){
       t.test(treated.cisplatin[x,], untreated.cisplatin[x,],paired= T)$p.value
   })
# 3. Wlicoxon rank sum test
pvalues.wilcoxon <- sapply(double.biomarker, function(x){
          wilcox.test(treated.cisplatin[x,], untreated.cisplatin[x,],paired= T)$p.value
      })
# Comparison t-tests
cbind(pvalues.welch, pvalues.wilcoxon)



#Step 2
#influence of cisplatin on the biomarkers gene expression in different cell lines
double.biomarker.FC = FC.cisplatin[double.biomarker,]
colfunc <- colorRampPalette(c("firebrick","firebrick3","lightcoral",
                              "lightyellow","lightskyblue1","steelblue1",
                              "steelblue3", "darkblue"))

#heatmap
pheatmap(double.biomarker.FC,
         color = colfunc(25),
         cluster_cols = TRUE,
         clustering_rows = TRUE,
         clustering_method ="ward.D2",
         treeheight_row = 20,
         treeheight_col = 30,
         cutree_cols = 5,
         cutree_rows = 2, 
         legend =TRUE,
         show_colnames = FALSE,
         legend_breaks = c(-3:3),
         legend_labels = c("red= downregulation","","","","","","blue= upregulation"),
         border_color = FALSE,
         scale = "column",
         main = "Influence of cisplatin on the biomarkers gene expression")

#heatmap with annotation
colnames(double.biomarker.FC) <- meta[95:149,2]
annotation = data.frame(Cancertype = cellline$Cancer_type)
rownames(annotation) = cellline$Cell_Line_Name


pheatmap(double.biomarker.FC,
         color = colfunc(25),
         cluster_cols = TRUE,
         clustering_rows = TRUE,
         clustering_method ="ward.D2",
         treeheight_row = 20,
         treeheight_col = 30,
         annotation_col = annotation,
         cutree_cols = 5,
         cutree_rows = 2, 
         legend =TRUE,
         show_colnames = FALSE, 
         legend_breaks = c(-3:3),
         legend_labels = c("red= downregulation","","","","","","blue= upregulation"),
         border_color = "white",
         scale = "column",
         main = "Influence of cisplatin on the biomarkers gene expression")


#checking if number of cluster is accurate (elbow plot) 
wss = sapply(2:8, function(k) {
  kmeans(x = t(double.biomarker.FC), centers = k)$tot.withinss
})
plot(2:8, wss, type = "b", pch = 19, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares")




##Step 3
#further analysis of biomarker
#amplification, deletation or neutral gene copy number
copynumber.biomarker = copynumber[double.biomarker,]
copynumber.biomarker = copynumber.biomarker[-11,]
copynumber.biomarker = as.matrix(copynumber.biomarker)


quantiles = as.matrix(quantile(copynumber.biomarker, probs= c(0, 0.2, 0.25, 0.5, 0.75, 0.8,1)))
plot(density(copynumber.biomarker))
abline(v= quantiles[c(2,6),1], col= c("red", "blue"), lty =2)
copynumber.quali = ifelse(copynumber.biomarker < (quantiles[2,1]), (-1), ifelse (copynumber.biomarker > (quantiles[6,1]), 1, 0))



#heatmap
colfunc2 <- colorRampPalette(c("firebrick2", "grey88", "deepskyblue3"))
colfunc2(3)

pheatmap(copynumber.quali,
         color = colfunc2(3),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "ward.D2",
         cutree_rows = 3,
         cutree_cols = 3, 
         legend = TRUE,
         legend_breaks = c(-3:3),
         legend_labels = c("red=deletion","","","", "", "", "blue=amplification"),
         scale = "column",
         border_color = "white",
         show_colnames = FALSE,
         main =  "Connection between biomarker and gene alterations")

#checking the optimal number of cluster (elbow plot) 
wss = sapply(2:7, function(k) {
  kmeans(x = t(copynumber.biomarker), centers = k)$tot.withinss
})
plot(2:7, wss, type = "b", pch = 19, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares")

kmeans= kmeans(x = t(copynumber.biomarker), centers = 3, nstart = 10)

#heatmap with annotation
pheatmap(copynumber.quali,
         color = colfunc2(3),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "ward.D2",
         cutree_cols = 3,
         cutree_rows = 3,
         annotation = annotation,
         legend = TRUE,
         legend_breaks = c(-3:3),
         legend_labels = c("red=deletion","","","", "", "", "blue=amplification"),
         scale = "column",
         border_color = "white",
         show_colnames = FALSE,
         main =  "Connection between biomarker and gene alterations")
