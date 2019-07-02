##loading data

library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)
results = readRDS(paste0(wd, "/data/results.RDS"))
results_cisplatin = subset (results , drug == "cisplatin")

# overview RGES values - damit wir besser einsch?tzen k?nnen was f?r uns gut oder schlecht ist
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
#nein es sind f?r beide drugs nur Ausrei?er - hat es dann vllt etwas mit den celllines zu tun? 
boxplot(results[which(results$cell == "UACC-62"),2], main = "RGES for cellline UACC-62 (Melanoma)",  ylim = c(min,max))
boxplot(results[which(results$cell == "OVCAR-4"),2], main = "RGES for cellline OVCAR-4 (Ovarian)",  ylim = c(min,max))
#auch hier sind es Ausrei?er. Das hei?t die guten RGES Werte liegen allein an der speziellen Kombination cellline+drug
#ich hab mal kurz gegoogelt und paclitaxel wird bei ovarienkarzinom verwendet, das macht also voll Sinn
#bortezomib wird nicht f?r melanome eingesetzt sonder f?r multiples myelom. da k?nnen wir dann sagen dass das ein Hinweis sein k?nnte
#das Medikament auch f?r Melanome zu testen


#Teresaaaa IC50 
library(reshape)

wd = dirname(rstudioapi::getSourceEditorContext()$path)
results = readRDS(paste0(wd, "/data/results.RDS"))
ic50 = readRDS(paste0(wd, "/data/NegLogGI50.RDS"))
meta = read.delim(paste0(wd, "/data/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t") 

#melt function
ic50 = t(ic50)
melt.data <- melt(ic50)
melt.data = as.matrix(melt.data)

#remove NAs
rmv.rows = apply(melt.data, 1, function(x) {
  sum(is.na(x))
})
melt.ic50 = melt.data[-which(rmv.rows > 0),]
rm(melt.data)

#celllines der Ic50 aussortieren - an RGES_results anpassen
colnames(melt.ic50) = c("cell", "drug", "IC50")

#zwei versuche die aber keinen Sinn machen
output.dataset = sapply(seq_along(melt.ic50), function(a) {
  out <- melt.ic50[,which(melt.ic50 == results$cell & melt.ic50 == results$drug)]
  return(out)
})

melt.ic50 = as.data.frame(melt.ic50)
output.dataset = sapply(1:894, function(a) {
  subset(melt.ic50, melt.ic50[a,1] == results$cell & melt.ic50$drug == results$drug)
})

##neuer Versuch Jojo
results.matrix.name = as.matrix(results)
melt.ic50.name = melt.ic50
rownames(results.matrix.name)= make.names(results.name[,4], unique = TRUE)
rownames(melt.ic50.name)=make.names(melt.ic50.name[,1], unique = TRUE)

names = rownames(results.matrix.name)
output.dataset <- sapply(seq_along(names), function(a) {
  name_picker <- names[a]
  out <- melt.ic50.name[,which(rownames(melt.ic50.name) == name_picker)]
  return(out)
})


########## drug efficacy plots 
#correlated to IC50

ic50 <- aggregate(standard_value ~ pert_iname, lincs_drug_activity_confirmed,median)

drug_activity_rges <- merge(results, ic50, by.x="cell", by.y="pert_iname")

drug_activity_rges <- aggregate(cbind(RGES, standard_value) ~ name, drug_activity_rges, median)

plot(drug_activity_rges$RGES, log(drug_activity_rges$standard_value, 10))
cor_test <- cor.test(drug_activity_rges$RGES, log(drug_activity_rges$standard_value, 10))

drug_activity_rges <- drug_activity_rges[order(drug_activity_rges$RGES),]

lm_cmap_ic50 <- lm(RGES ~ log(standard_value, 10), drug_activity_rges)


pdf(paste( "fig/", cancer, "rges_ic50_cmap_data_", landmark, ".pdf", sep=""))
ggplot(drug_activity_rges, aes(RGES, log(drug_activity_rges$standard_value, 10)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  annotate("text", label = paste(cancer, ",", "MCF7", sep=""), 
           x = 0, y = 8.1, size = 6, colour = "black") +
  annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=2), sep=""), 
           x = 0, y = 7.7, size = 6, colour = "black") +
  annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 7.3, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab("RGES") + guides(shape=FALSE, size=FALSE) +
  ylab("log10(IC50) nm") + coord_cartesian(xlim = c(-0.5, 0.5), ylim=c(-1, 8)) 
dev.off()


