### IC50 spalte an results anf?gen -> drug_activity_rges
##loading data
library(reshape)

wd = dirname(rstudioapi::getSourceEditorContext()$path)
results = readRDS(paste0(wd, "/data/results.RDS"))
ic50 = readRDS(paste0(wd, "/data/NegLogGI50.RDS"))
meta = read.delim(paste0(wd, "/data/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t") 

##melt function
ic50 = t(ic50)
melt.data <- melt(ic50)
melt.data = as.matrix(melt.data)
melt.ic50 = as.data.frame(melt.data)
colnames(melt.ic50) = c("cell", "drug", "IC50")

##remove NAs
rmv.rows = apply(melt.data, 1, function(x) {
  sum(is.na(x))
})
which(rmv.rows > 0)
melt.ic50 = melt.data[-which(rmv.rows > 0),]
rm(melt.data)

##Cellines der IC50 aussortieren um sie an RGES_results anpassen
IC50.value = c(rep(as.numeric(0),819))
results = cbind(results, IC50.value)

i=1
j=1

while(i<895)
{while(j<820)
{
  if(isTRUE(melt.ic50[i,1]== results[j,4])
     & (melt.ic50[i,2] == results[j,5]))
  {results[j,9] = as.numeric(melt.ic50[i,3])
  }
  j = j +1
}
  j= 1
  i=i+1}

which(results$IC50.value == 0)
drug_activity_rges = results[-which(results$IC50.value == 0),]


saveRDS(drug_activity_rges, file = "drug_activity_rges.rds")

# limit data set to numerical values
rges.ic50 <- as.data.frame(drug_activity_rges[,c(2,9)])

# compute correlation
cor(rges.ic50$RGES, rges.ic50$IC50.value, method = "spearman")
# check significance
cor.test(rges.ic50$RGES, rges.ic50$IC50.value, method = "spearman")
#plot correlation
library(ggplot2)
ggplot(rges.ic50, aes(rges.ic50$RGES, rges.ic50$IC50.value)) +
  geom_point(color = "blue", size = 1) +
  scale_size(range = c(2,5)) +
  xlab("RGES") + 
  ylab("IC50 [nm]") +
  xlim(0.1,0.25)+
  ylim(-1e+09,0.2e+09)

# univariate regression model
# split data set in test and training set
# take 200 random samples to form the training set
i.train = sample(1:nrow(rges.ic50), 730)
rges.ic50.train = rges.ic50[i.train, ]
rges.ic50.test = rges.ic50[-i.train, ]
#learn lm on training set
lm.rges_ic50 = lm(IC50.value ~ RGES, data = rges.ic50.train)
summary(lm.rges_ic50)

# normal distribution of residuals?
plot(lm.rges_ic50, which = c(1), pch=20, col="blue", main = "Scatter plot Residuals - Fitted values")
plot(lm.rges_ic50, which = c(2), pch=20, col="blue", main = "QQ plot Residuals")
## correlation residuals x-values?
cor(rges.ic50.train$RGES, lm.rges_ic50$residuals)

# use model to predict ic50 by rges
pred = predict(lm.rges_ic50, newdata = rges.ic50.test)
plot(rges.ic50.test$IC50.value, pred, xlab = "Real Values", ylab = "Predicted Values", pch=20, col="blue")
abline(0, 1, col = "red")
#compute RMSE, check validility
n = nrow(rges.ic50.train)
rmse.train = sqrt(1/n * sum(lm.rges_ic50$residuals^2))
n = nrow(rges.ic50.test)
residuals = rges.ic50.test$IC50.value - pred
rmse.test = sqrt(1/n * sum(residuals^2))


#Multiple Regression
# 1. Include biomarkers to RGES matrix only for cisplatin!
double.biomarker.FC = readRDS(paste0(wd, "/data/double.biomarker.FC.RDS"))
drug_activity_rges.cisplatin = subset (drug_activity_rges , drug == "cisplatin")
# samples as rows
biomarker.FC = t(double.biomarker.FC)

# fit both matrices
biomarker.FC.fit = subset(biomarker.FC, rownames(biomarker.FC) %in% drug_activity_rges.cisplatin$sample)

# limit the data set to numerical values: rges, ic50 and biomarkers:
rges.ic50.biomarkers <- as.data.frame(cbind(drug_activity_rges.cisplatin[,c(2,9)], biomarker.FC.fit))


# correlation between biomarkers, rges and ic50 is visualized
# produce pairwise scatter plots
pairs(rges.ic50.biomarkers, col = "blue", pch = 20, main = "Scatterplots RGES, IC50, Biomarkers")
## matrix of correlations
cor.mat = cor(rges.ic50.biomarkers, method = "spearman")
heatmap(cor.mat, col = cm.colors(256), symm = T, main = "Correlation Heatmap RGES, IC50 and Biomarkers")

# multiple regression model
# create training and test set
train.multiple = sample(1:nrow(rges.ic50.biomarkers), 45)
train.set.multiple = rges.ic50.biomarkers[train.multiple, ]
test.set.multiple = rges.ic50.biomarkers[-train.multiple, ]

# learn model, all data included
model.multiple = lm(IC50.value ~ ., data = train.set.multiple)
summary(model.multiple)

# prove residuals
plot(model.multiple, which = c(1), col="blue", pch = 20, main = "Scatterplot Residuals - Fitted values")
plot(model.multiple, which = c(2), col="blue", pch = 20, main = "QQ plot Residuals")
cor(rges.ic50.train[,-2], lm.rges_ic50$residuals)

# predict ic50
predict.multiple = predict(model.multiple, newdata = test.set.multiple)
plot(test.set.multiple$IC50.value, predict.multiple, xlab = "Real Values", ylab = "Predicted Values", pch=20, col="blue")
abline(0, 1, col = "red")

#computation of RMSE 
n = nrow(train.set.multiple)
rmse.train = sqrt(1/n * sum(model.multiple$residuals^2))
n = nrow(test.set.multiple)
residuals = test.set.multiple$IC50.value - predict.multiple
rmse.test = sqrt(1/n * sum(residuals^2))
rmse.train
rmse.test

# 2. model with specific biomarkers (PTPRG, COMMD10, GMDS, LRBA)
#same training and test set is used
model.multiple.biomarkers = lm(IC50.value ~ PTPRG + COMMD10 + GMDS + LRBA, data = train.set.multiple)
summary(model.multiple.biomarkers)

# prove residuals
plot(model.multiple.biomarkers, which = c(1), col = "blue", pch = 20, main = "Scatterplot Residuals - Fitted values")
plot(model.multiple.biomarkers, which = c(2), col = "blue", pch = 20, main = "QQ plot Residuals")

# predict ic50
predict.biomarkers = predict(model.multiple.biomarkers, newdata = test.set.multiple)
plot(test.set.multiple$IC50.value, predict.biomarkers, xlab = "Real Values", ylab = "Predicted Values", pch=20, col="blue")
abline(0, 1, col = "red")

# compute RMSE
n = nrow(train.set.multiple)
rmse.train = sqrt(1/n * sum(model.multiple.biomarkers$residuals^2))
n = nrow(test.set.multiple)
residuals = test.set.multiple$IC50.value - predict.biomarkers
rmse.test = sqrt(1/n * sum(residuals^2))

# for a better model, pca can be used
pca = prcomp(rges.ic50.biomarkers[, -2])
barplot(pca$rotation[, 1], horiz = TRUE, main = "PC1", col = "lightblue", las=1)

# compute model with pcas instead of original variables
model.pca = lm(rges.ic50.biomarkers$IC50.value ~ pca$x)
summary(model.pca)
plot(model.pca, which = c(1), col = "blue", pch = 20, main = "Scatterplot Residuals - Fitted values")
plot(model.pca, which = c(2), col = "blue", pch = 20, main = "QQ plot Residuals")
# do the PCs correlate?
cor.pca = cor(pca$x)
heatmap(cor.pca, col = cm.colors(256), main = "Heatmap Correlation PCs")
# would it be better to redo the model with only a few of PCs?
plot(pca, type ="lines", main = "Elbow plot of PCA")
