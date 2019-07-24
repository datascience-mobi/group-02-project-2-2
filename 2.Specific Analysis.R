## packages 
library(rstudioapi)
library(gplots)
library(pheatmap)
library(viridis)
library(VennDiagram)

## necessary data from the general analysis
wd = dirname(rstudioapi::getSourceEditorContext()$path)
meta = read.delim(paste0(wd, "/data/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t") 
untreated = readRDS(paste0(wd, "/data/NCI_TPW_gep_untreated.RDS"))
treated = readRDS(paste0(wd, "/data/NCI_TPW_gep_treated.RDS"))
copynumber = readRDS(paste0(wd, "/data/CCLE_copynumber.RDS"))
cellline = read.delim(paste0(wd, "/data/cellline_annotation.tsv"), header = TRUE, sep = "\t")

#scale
treated.scaled <- scale(treated)
untreated.scaled <- scale(untreated)

#FC
log2FC.treated.untreated <- (treated.scaled-untreated.scaled) #values are already log2 transformed

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

mean.FC.genes <- apply(FC.cisplatin, 1, mean)
mean.FC.celllines <- apply(FC.cisplatin, 2, mean)

#removing unnecessary data, functions...
rm(cisplatin.col, j, treated, untreated, log2FC.treated.untreated)

#standard deviation cisplatin
sd.cisplatin = apply(FC.cisplatin, 1, sd)
mean.FC.genes.sd = as.data.frame(cbind(mean.FC.genes, sd.cisplatin))
mean.FC.ordered = mean.FC.genes.sd[order(mean.FC.genes.sd$mean.FC.genes, decreasing=TRUE), ]
plot(density(sd.cisplatin))

#finding/visualizing most extreme FC values for cisplatin
lowest.FC = mean.FC.ordered[13280:13299,]
par(mar = c(5, 7, 5, 3))
barplot(lowest.FC$mean.FC.genes,
        horiz = TRUE,
        xlim = c(-1.0, 0),
        main= "lowest log2 FC-values for cisplatin",
        xlab= "mean log2FC values in different celllines",
        names.arg = rownames(lowest.FC),
        col= "firebrick",
        las=1,
        border = "white", 
        cex.names =0.8)


highest.FC = mean.FC.ordered[1:20,]
par(mar = c(5, 10, 5, 3))
barplot(highest.FC$mean.FC.genes,
        horiz = TRUE,
        xlim = c(0, 1),
        main= "highest log2 FC-values for cisplatin",
        xlab= "mean log2FC values in different celllines",
        names.arg = rownames(highest.FC),
        col= "lightgreen",
        las=1,
        border = "white",
        cex.names =0.8)


#matrix containing the biomarker found through the FC 
highest.FC = as.matrix(highest.FC) 
lowest.FC = as.matrix(lowest.FC)
biomarker1.FC = as.matrix(rbind(highest.FC, lowest.FC))
highest.names <- row.names(highest.FC)
lowest.names <- row.names(lowest.FC)
row.names(biomarker1.FC) <- c(highest.names, lowest.names)
biomarker1 <- c(highest.names, lowest.names)

#sort out biomarker that don?t change in the same "direction" for most cell lines
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

rm(is.neg, a, i, j, biomarker2.down, biomarker2.up, highest.names, lowest.names)

grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1 = 40,
  area2 = 668,
  cross.area = 15,
  fill = c("yellowgreen", "steelblue2"),
  category = c("FC-Criterium", "Consistent up- or downregulation"),
  lty = "blank",
  cex = 2,
  cat.cex = 1.5,
  cat.pos = c(0, 3),
  cat.dist = -0.5,
  cat.just = list(c(-0.0, 11), c(1,1.5)),
  cat.col = c("yellowgreen", "steelblue2"), 
  ext.pos = c(-3,2),
  ext.dist = c(-0.4, -0.05),
  ext.length = 0.8
)
grid.draw(venn.plot)

rm(venn.plot)
## Step 1 (b)
# ttest to verify significance of the biomarker
#Creating seperated matrices containing cisplatin
treated.cisplatin <- treated.scaled[,grep ("cisplatin", colnames(treated.scaled))]
untreated.cisplatin <- untreated.scaled[,grep("cisplatin", colnames(treated.scaled))]

# 1. Normality was checked for all biomarkers through QQplots. This is an example.
par(mar = c(5, 5, 5, 5))
qqnorm(treated.cisplatin["POLR3B", ], main = "QQplot to check normality of the POLR3B gene")
qqline(treated.cisplatin["POLR3B", ])



# 2. Welch two sample t-test. Check if the expressional change of biomarkers is significant. 
pvalues.welch <- sapply(double.biomarker, function(x){
       t.test(treated.cisplatin[x,], untreated.cisplatin[x,],paired= T)$p.value
   })
plot(density(pvalues.welch),  main = "Welch two sample t-test - density of p values")

# Standard deviation
sd.cisplatin = apply(FC.cisplatin, 1, sd)
mean.FC.genes.sd = as.data.frame(cbind(mean.FC.genes, sd.cisplatin))
mean.FC.ordered = mean.FC.genes.sd[order(mean.FC.genes.sd$mean.FC.genes, decreasing=TRUE), ]
head(mean.FC.genes.sd)
plot(density(sd.cisplatin))


#checking consistency with standard deviation for double.biomarker
double.biomarker.FC = FC.cisplatin[double.biomarker,]
mean.double <- apply(double.biomarker.FC, 1, mean)
sd.double <- apply(double.biomarker.FC, 1, sd)
double = cbind(mean.double, sd.double)
colnames(double) = c("mean FC","SD" )

#Step 2
#influence of cisplatin on the biomarkers gene expression in different cell lines
double.biomarker.FC = FC.cisplatin[double.biomarker,]
colfunc <- colorRampPalette(c("firebrick","firebrick3","lightcoral","lightyellow","lightskyblue1","steelblue1","steelblue3", "darkblue"))


#checking how many clusters we will need in the heatmap (elbow plot) 
wss = sapply(2:8, function(k) {
  kmeans(x = t(double.biomarker.FC), centers = k)$tot.withinss
})
plot(2:8, wss, type = "b",
     pch = 19,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares",
     main = "Elbow plot for kmeans clustering of celllines")

#heatmap with annotation
colnames(double.biomarker.FC) <- meta[95:149,2]
annotation = data.frame(Cancertype = cellline$Cancer_type)
rownames(annotation) = cellline$Cell_Line_Name

par(mar = c(10,20,10,5))
pheatmap(double.biomarker.FC,
         color = colfunc(25),
         cluster_cols = TRUE,
         clustering_rows = TRUE,
         clustering_method ="ward.D2",
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_col = annotation,
         legend = T,
         legend_breaks = c(-1,1),
         legend_labels = c("up", "down"),
         show_colnames = F,
         cutree_rows = 2,
         cutree_cols = 4,
         border_color = "white",
         scale = "column",
         main = "Influence of cisplatin on the biomarkers gene expression")



##Step 3
#further analysis of biomarker
#amplification, deletation or neutral gene copy number
copynumber.biomarker = as.matrix(copynumber[double.biomarker,])
copynumber.quali = ifelse(copynumber.biomarker <= (-1), (-1), ifelse (copynumber.biomarker >= (1), 1, 0))

#distribution of genetic alterations
d <- density(copynumber.biomarker)
plot(d, main = "Distribution of genetic alterations in biomarker genes")
abline(v= c(-1,1), col= c("red", "blue"), lty =2)
legend(-6, 1.2, legend = c("Deletion", "Amplification"), col = c("red", "blue"), lty = 2)

#checking the optimal number of cluster (elbow plot) 
wss = sapply(2:7, function(k) {
  kmeans(x = t(copynumber.quali), centers = k)$tot.withinss
})
plot(2:7, wss, 
     type = "b", pch = 19, 
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares",
     main = "Elbowplot kmeans clustering - biomarker gene alterations")

kmeans= kmeans(x = t(copynumber.quali), centers = 4, nstart = 10)
kmeans$totwithinss

#heatmap
colfunc2 <- colorRampPalette(c("firebrick2", "grey88", "deepskyblue3"))

colnames(double.biomarker.FC) <- meta[95:149,2]
annotation = data.frame(Cancertype = cellline$Cancer_type)
rownames(annotation) = cellline$Cell_Line_Name

pheatmap(copynumber.quali,
         color = colfunc2(3),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "ward.D2",
         legend = TRUE,
         legend_breaks = c(-1,0,1),
         legend_labels = c("deletion","", "amplification"),
         border_color = "white",
         show_colnames = FALSE,
         cutree_rows = 4,
         cutree_cols = 3,
         treeheight_row = 20,
         treeheight_col = 20,
         annotation = annotation)