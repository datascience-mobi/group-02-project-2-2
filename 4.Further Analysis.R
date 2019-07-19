##loading data

library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)
results = readRDS(paste0(wd, "/data/results.RDS"))
results_cisplatin = subset (results , drug == "cisplatin")

###RGES analysis
## overview RGES values
quantile(results$RGES)
quantile(results_cisplatin$RGES)
min = min(results$RGES)
max = max(results$RGES)

#plot
quantile.res = as.matrix(quantile(results$RGES))
plot(density(results$RGES))
abline(v= quantile.res[c(2,4),1], col= c("red", "blue"), lty =2)
legend(-6, 1.2, legend = c("25 % quantile", "75% quantile"), col = c("red", "blue"), lty = 1:2)

quantile.cis = as.matrix(quantile(results_cisplatin$RGES))
plot(density(results_cisplatin$RGES))
abline(v= quantile.cis[c(2,4),1], col= c("red", "blue"), lty =2)
legend(-6, 1.2, legend = c("25 % quantile", "75% quantile"), col = c("red", "blue"), lty = 1:2)


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

# check for negative RGES values again
results_neg =results[which(results$RGES < 0),]

#distribution RGES in different drugs
dev.new(width=5, height=4)
par(mfrow=c(5,3), mar=c(2,2,2,2))
sapply(1:length(drug),
       function(x) boxplot(results[which(results$drug == drug[x]),2],
                           main = drug[x], 
                           ylim = c(min,max)))