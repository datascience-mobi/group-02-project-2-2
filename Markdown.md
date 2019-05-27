# project-02-group-02
group-02-project-2-2 created by GitHub Classroom 
# Chunk options
```{r setup, include=FALSE}
knitr::opts_chunks$set(echo = TRUE)
knitr::opts_chunk$set(cache = TRUE)

#Load Data  
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)

#rename 
```{r read_data}
meta = read.delim(paste0(wd, "/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t")  
untreated = readRDS(paste0(wd, "/NCI_TPW_gep_untreated.RDS"))  
treated = readRDS(paste0(wd, "/NCI_TPW_gep_treated.RDS"))  
ic50 = readRDS(paste0(wd, "/NegLogGI50.RDS"))  
basalexp = readRDS(paste0(wd, "/CCLE_basalexpression.RDS"))  
copynumber = readRDS(paste0(wd, "/CCLE_copynumber.RDS"))  
mutations = readRDS(paste0(wd, "/CCLE_mutations.RDS"))  
cellline = read.delim(paste0(wd, "/cellline_annotation.tsv"), header = TRUE, sep = "\t")  
drug = read.delim(paste0(wd, "/drug_annotation.tsv"), header = TRUE, sep = "\t")  
```

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


```{r}
library(rstudioapi) 
wd = dirname(rstudioapi::getSourceEditorContext()$path)
```
Now we will rename our data to get a better overview.
```{r read_data}
meta = read.delim(paste0(wd, "/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t") 
untreated = readRDS(paste0(wd, "/NCI_TPW_gep_untreated.RDS"))
treated = readRDS(paste0(wd, "/NCI_TPW_gep_untreated.RDS"))
ic50 = readRDS(paste0(wd, "/NegLogGI50.RDS"))
basalexp = readRDS(paste0(wd, "/CCLE_basalexpression.RDS"))
copynumber = readRDS(paste0(wd, "/CCLE_copynumber.RDS"))
mutations = readRDS(paste0(wd, "/CCLE_mutations.RDS"))
cellline = read.delim(paste0(wd, "/cellline_annotation.tsv"), header = TRUE, sep = "\t")
drug = read.delim(paste0(wd, "/drug_annotation.tsv"), header = TRUE, sep = "\t")
```
Let's check how many ``` NA ``` values we have.
```{r}
