##loading data

library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)
results = readRDS(paste0(wd, "/data/results.RDS"))

# overview RGES values - damit wir besser einschätzen können was für uns gut oder schlecht ist
quantile(results$RGES)
quantile(results_cisplatin$RGES)

# RGES cisplatin in general
results_cisplatin = subset (results , drug == "cisplatin")
boxplot(results_cisplatin$RGES, main = "RGES for cisplatin in different celllines")

# RGES cisplatin in different tissues 
tissue = c( "Renal", "Lung" , "Breast" , "Colon", "Prostate" , "Leukemia", "Ovarian", "Melanoma", "CNS")
mean_tissue =sapply(1:length(tissue), function(x) mean(subset(results_cisplatin$RGES , tissue == tissue[x])))
barplot(mean_tissue , name = tissue , las = 2 , horiz = FALSE ,col= "firebrick", border = "white" , main = "RGES for cisplatin in different tissues" , ylab = "RGES")

# find negative RGES values
results_neg = results[which(results$RGES < 0),]




