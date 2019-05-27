---
output: html_document
editor_options: 
  chunk_output_type: inline
---
# project-02-group-02
group-02-project-2-2 created by GitHub Classroom  

# Load Data  
```
library(rstudioapi)
wd = dirname(rstudioapi::getSourceEditorContext()$path)
```

## rename 
```{r}
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
## general overview
``` {r}
boxplot(treated, xlab = "treated", horizontal = F, border=cb)
boxplot(treated.scaled, xlab="treated scaled", horizontal =F, border=cb)
```