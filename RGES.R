### package for nbinomTest
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq")

###data
library(rstudioapi)
wd = dirname(rstudioapi::getSourceEditorContext()$path)
meta = read.delim(paste0(wd, "/data/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t") 
untreated = readRDS(paste0(wd, "/data/NCI_TPW_gep_untreated.RDS"))
treated = readRDS(paste0(wd, "/data/NCI_TPW_gep_treated.RDS"))
basalexp = readRDS(paste0(wd, "/data/CCLE_basalexpression.RDS"))
library(DESeq)


untreated.fit = subset(untreated, rownames(untreated) %in% rownames(basalexp))
treated.fit = subset(treated, rownames(treated) %in% rownames(basalexp))

#drug signature DSEq
mode(treated.fit) <- "integer"
mode(untreated.fit) <- "integer"
treated.untreated <- cbind(treated.fit,untreated.fit)
cds = newCountDataSet(countData = treated.untreated, conditions = c(rep("treated",819),rep("untreated",819)))
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
str( fitInfo(cds) )
plotDispEsts(cds)
dispersion.values = fData(cds)
nbinom.treated.untreated = nbinomTest(cds, "treated", "untreated")
plotMA(nbinom.treated.untreated)
hist(nbinom.treated.untreated$pval, breaks=100, col="skyblue", border="slateblue", main="p-values nbinom")



###basal anpassen an untreated - die neue heisst jetzt basal.fitted.untreated
# nurnoch die Gene dalassen welche in basal und untreated sind
basal.fit = subset(basalexp, rownames(basalexp) %in% rownames(untreated))

meta.matrix <- as.matrix(meta)
dataset.basal <- as.matrix(basal.fit)

new.basal.names <- as.character(meta.matrix[1:819,2])
output.dataset <- sapply(seq_along(new.basal.names), function(a) {
  name_picker <- new.basal.names[a]
  out <- dataset.basal[,which(colnames(dataset.basal) == name_picker)]
  return(out)
})

#als Matrix umformatieren und umnennen damit wir es erkennen
basal.fitted.untreated <- matrix(unlist(output.dataset), nrow = 11461, ncol = 819, byrow=FALSE, dimnames = NULL)
colnames(basal.fitted.untreated) <- make.names(new.basal.names, unique = TRUE)

# rownames: gene einfuegen
rownames(basal.fitted.untreated)= make.names(rownames(basal.fit), unique = TRUE)



#disease signature
#gleiches prinzip wie oben, jetzt zwischen basal.fitted.untreated und untreated.fit
mode(basal.fitted.untreated) <- "integer"
mode(untreated.fit) <- "integer"
basal.untreated <- cbind(basal.fitted.untreated, untreated.fit)
cds = newCountDataSet(countData = basal.untreated, conditions = c(rep("basal",819),rep("untreated",819)))
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
str( fitInfo(cds) )
plotDispEsts(cds)
dispersion.values2 <- fData(cds)
nbinom.basal.untreated = nbinomTest(cds, "basal", "untreated")
plotMA(nbinom.basal.untreated)
hist(nbinom.basal.untreated$pval, breaks=100, col="skyblue", border="slateblue", main="p-values nbinom")



## signature genes nach den Kriterien von Bin Chen
  #for the RGESexample code they use 978 genes, so we aim for the same number of genes
  #log2 Kriterium ganz raus weil dafuer sind unsere Werte viel zu klein

#disease signature
dz_signature <- subset(nbinom.treated.untreated, !is.na(padj) & !is.na(id) & id !='?' & padj < 0.5  & abs(log2FoldChange) != Inf )
dim(dz_signature)
gene.list.1 = c(dz_signature[,1])

#drugsignature
dr_signature <- subset(nbinom.basal.untreated, !is.na(padj) & !is.na(id) & id !='?' & padj < 0.5  & abs(log2FoldChange) != Inf )
dim(dr_signature)
gene.list.2 = c(dr_signature[,1])

#Liste von Genen, die nach Anwendung der Kriterien noch da sind 
gene.list.final = Reduce(intersect, list(gene.list.1, gene.list.2))
length(gene.list.final)


##log2 Fold-change and scaling 
untreated.fit = subset(untreated, rownames(untreated) %in% rownames(basalexp))
treated.fit = subset(treated, rownames(treated) %in% rownames(basalexp))

treated.fit <- scale(treated.fit)
untreated.fit <- scale(untreated.fit)
basal.fitted.untreated <- scale(basal.fitted.untreated)

drug_signature = log2(treated.fit/untreated.fit)
is.nan.data.frame <- function(x)     
  do.call(cbind, lapply(x, is.nan))
drug_signature[is.nan(drug_signature)] <- 0
dim(drug_signature)

disease_signature <- log2(untreated.fit/basal.fitted.untreated)
is.nan.data.frame <- function(x)     
  do.call(cbind, lapply(x, is.nan))
disease_signature[is.nan(disease_signature)] <- 0
dim(drug_signature)


#only keep FC values for dz_signature and dr_signature genes (also die die die Kriterien erfüllen )
drug_signature = subset(drug_signature, rownames(drug_signature) %in% gene.list.final)
dim(drug_signature)
disease_signature = subset(disease_signature, rownames(disease_signature) %in% gene.list.final)
dim(disease_signature)

## function from Bin Chen to compute RGES
cmap_score_new <- function(sig_up, sig_down, drug_signature) {
  num_genes <- nrow(drug_signature)
  ks_up <- 0
  ks_down <- 0
  connectivity_score <- 0
  drug_signature[,"rank"] <- rank(drug_signature[,"rank"])
  up_tags_rank <- merge(drug_signature, sig_up, by.x = "ids", by.y = 1)
  down_tags_rank <- merge(drug_signature, sig_down, by.x = "ids", by.y = 1)
  up_tags_position <- sort(up_tags_rank$rank)
  down_tags_position <- sort(down_tags_rank$rank)
  num_tags_up <- length(up_tags_position)
  num_tags_down <- length(down_tags_position)
  if(num_tags_up > 1) {
    a_up <- 0
    b_up <- 0
    a_up <- max(sapply(1:num_tags_up,function(j) {
      j/num_tags_up - up_tags_position[j]/num_genes
    }))
    b_up <- max(sapply(1:num_tags_up,function(j) {
      up_tags_position[j]/num_genes - (j-1)/num_tags_up
    }))
    if(a_up > b_up) {
      ks_up <- a_up
    } else {
      ks_up <- -b_up
    }
  }else{
    ks_up <- 0
  }
  if (num_tags_down > 1){
    a_down <- 0
    b_down <- 0
    a_down <- max(sapply(1:num_tags_down,function(j) {
      j/num_tags_down - down_tags_position[j]/num_genes
    }))
    b_down <- max(sapply(1:num_tags_down,function(j) {
      down_tags_position[j]/num_genes - (j-1)/num_tags_down
    }))
    if(a_down > b_down) {
      ks_down <- a_down
    } else {
      ks_down <- -b_down
    }
  }else{
    ks_down <- 0
  }
  if (ks_up == 0 & ks_down != 0){ #only down gene inputed
    connectivity_score <- -ks_down 
  }else if (ks_up !=0 & ks_down == 0){ #only up gene inputed
    connectivity_score <- ks_up
  }else if (sum(sign(c(ks_down,ks_up))) == 0) {
    connectivity_score <- ks_up - ks_down # different signs
  }else{
    connectivity_score <- ks_up - ks_down
  }
  return(connectivity_score)
}



##CLAUDIA 
sig.ids <- c(1:819) ### cell lines bekommen sig.ID: sig.ID=nrow(meta[,1:819])=ncol(treated usw)
dz_genes_up <- subset(nbinom.basal.untreated, log2FoldChange > 0,select="id")
##bei uns ist die ID der Name des Gens (nicht verwechseln mit sig.ID)
dz_genes_down <- subset(nbinom.basal.untreated, log2FoldChange <0 ,select="id")
gene.list <- rownames(disease_signature)

dz_cmap_scores <- NULL
count <- 0
for(count in sig.ids){
  print(count)
  cmap_exp_signature <- data.frame(gene.list,  rank(-1 * disease_signature[, count], ties.method="random"))    
  colnames(cmap_exp_signature) <- c("ids","rank") 
  dz_cmap_scores <- c(dz_cmap_scores, cmap_score_new(dz_genes_up,dz_genes_down,cmap_exp_signature)) #hier wird der RGES berechnet
  count <- count + 1
}

results <- data.frame(sig.ID = sig.ids, RGES = dz_cmap_scores)
results <- cbind(results, meta[1:819,]) 

