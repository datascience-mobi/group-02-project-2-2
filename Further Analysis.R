##loading data

library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)
results = readRDS(paste0(wd, "/data/results.RDS"))
results_cisplatin = subset (results , drug == "cisplatin")

# overview RGES values - damit wir besser einschätzen können was für uns gut oder schlecht ist
quantile(results$RGES)
quantile(results_cisplatin$RGES)
min = min(results$RGES)
max = max(results$RGES)

# RGES cisplatin in general
boxplot(results_cisplatin$RGES, main = "RGES for cisplatin in different celllines")

# mean RGES cisplatin in different tissues 
tissue = c( "Renal", "Lung" , "Breast" , "Colon", "Prostate" , "Leukemia", "Ovarian", "Melanoma", "CNS")
mean_tissue =sapply(1:length(tissue), function(x) mean(subset(results_cisplatin$RGES , tissue == tissue[x])))
barplot(mean_tissue , name = tissue , las = 2 , horiz = FALSE ,col= "firebrick", border = "white" , main = "mean RGES for cisplatin in different tissues" , ylab = "RGES")
nrow(min(results_cisplatin$RGES))


# mean RGES for different drugs
drug = c("5-Azacytidine", "bortezomib", "cisplatin","dasatinib","doxorubicin","erlotinib","geldanamycin","gemcitibine","lapatinib","paclitaxel","sirolimus","sorafenib","sunitinib","topotecan","vorinostat")
mean_drug =sapply(1:length(drug), function(x) mean(subset(results$RGES , drug == drug[x])))
barplot(mean_drug , name = drug , las = 2 , horiz = FALSE ,col= "forestgreen", border = "white" , main = "mean RGES for different drugs" , ylab = "RGES")

#distribution RGES in different drugs
#macht ein neues Fenster auch (dev.new) weil ich finde man erkennt das in dem plot Fenster nicht
dev.new(width=5, height=4)
par(mfrow=c(5,3), mar=c(2,2,2,2))
sapply(1:length(drug),
       function(x) boxplot(results[which(results$drug == drug[x]),2],
                           main = drug[x], 
                           ylim = c(min,max)))

# find negative RGES values
results_neg =results[which(results$RGES < 0),]
#sind diese drugs in allen celllines gut?
boxplot(results[which(results$drug == "bortezomib"),2], main = "RGES for bortezomib in different celllines",  ylim = c(min,max))
boxplot(results[which(results$drug == "paclitaxel"),2], main = "RGES for paclitaxel in different celllines",  ylim = c(min,max))
#nein es sind für beide drugs nur Ausreißer - hat es dann vllt etwas mit den celllines zu tun? 
boxplot(results[which(results$cell == "UACC-62"),2], main = "RGES for cellline UACC-62 (Melanoma)",  ylim = c(min,max))
boxplot(results[which(results$cell == "OVCAR-4"),2], main = "RGES for cellline OVCAR-4 (Ovarian)",  ylim = c(min,max))
#auch hier sind es Ausreißer. Das heißt die guten RGES Werte liegen allein an der speziellen Kombination cellline+drug
#ich hab mal kurz gegoogelt und paclitaxel wird bei ovarienkarzinom verwendet, das macht also voll Sinn
#bortezomib wird nicht für melanome eingesetzt sonder für multiples myelom. da können wir dann sagen dass das ein Hinweis sein könnte
#das Medikament auch für Melanome zu testen



