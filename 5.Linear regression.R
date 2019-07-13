### load the data
wd = dirname(rstudioapi::getSourceEditorContext()$path)
drug_activity_rges = readRDS(paste0(wd, "/data/drug_activity_rges.RDS"))


# limit the data set to numerical values, convert IC50 values as they are in -log10 scale
#keine ahnung ob das sinn macht, was ich hier getan habe?
ic50 = drug_activity_rges[,9]
IC50.value = 10^-(ic50)
RGES = drug_activity_rges[,2]
rges.ic50 <- as.data.frame(cbind(RGES, IC50.value))

# compute correlation
cor(rges.ic50$RGES, rges.ic50$IC50.value)
# check significance
cor.test(rges.ic50$RGES, rges.ic50$IC50.value)

## take 200 random samples to form the training set
i.train = sample(1:nrow(rges.ic50), 730)
## 
rges.ic50.train = rges.ic50[i.train, ]
rges.ic50.test = rges.ic50[-i.train, ]
#learn lm on training set
l.train = lm(IC50.value ~ RGES, data = rges.ic50.train)
summary(l.train)
plot(l.train)

# normal distribution of residuals?
hist(l.train$residuals, breaks = 20)
qqnorm(l.train$residuals)
qqline(l.train$residuals)
## correlation residuals x-values?
cor(rges.ic50.train$RGES, l.train$residuals)

# use model to predict ic50 by rges
pred = predict(l.train, newdata = rges.ic50.test)

#compute RMSE
n = nrow(rges.ic50.train)
rmse.train = sqrt(1/n * sum(l.train$residuals^2))
rmse.train
n = nrow(rges.ic50.test)
residuals = rges.ic50.test$IC50.value - pred
rmse.test = sqrt(1/n * sum(residuals^2))
rmse.test




#Multiple Regression
# 1. Include biomarkers to RGES matrix only for cisplatin!
double.biomarker.FC = readRDS(paste0(wd, "/data/double.biomarker.FC.RDS"))

drug_activity_rges.cisplatin = subset (drug_activity_rges , drug == "cisplatin")
# transformieren, damit samples in zeile
biomarker.FC = t(double.biomarker.FC)

# Problem: dimensionen der matrizen unterschiedlich, drugactivity eine zeile weniger als biomarker.FC scheinbar in Zeile 26 unterschiedlich?
# Anpassen der beiden matrizen aneindander:
biomarker.FC.fit = subset(biomarker.FC, rownames(biomarker.FC) %in% drug_activity_rges.cisplatin$sample)
results_biomarkers = cbind(drug_activity_rges.cisplatin, biomarker.FC.fit)

# limit the data set to numerical values: rges, ic50 and biomarkers:
# ic50 values are in -log10 scale
ic50.cisplatin = drug_activity_rges.cisplatin[,9]
IC50.value.cisplatin = 10^-(ic50.cisplatin)
RGES.cisplatin = drug_activity_rges.cisplatin[,2]
rges.ic50.biomarkers <- as.data.frame(cbind(RGES.cisplatin, IC50.value.cisplatin,biomarker.FC.fit))

# correlation between biomarkers, rges and ic50 is visualized
# produce pairwise scatter plots
pairs(rges.ic50.biomarkers, col = "blue", pch = 20)
## matrix of correlations
cor = cor(rges.ic50.biomarkers)
heatmap(cor, col = cm.colors(256))

# multiple regression model
# create training and test set
train.multiple = sample(1:nrow(rges.ic50.biomarkers), 45)
## 
train.set.multiple = rges.ic50.biomarkers[train.multiple, ]
test.set.multiple = rges.ic50.biomarkers[-train.multiple, ]

# learn model, all data included
model.multiple = lm(IC50.value.cisplatin ~ ., data = train.set.multiple)
summary(model.multiple)

plot(model.multiple)

predict.multiple = predict(model.multiple, newdata = test.set.multiple)

#computation of RMSE necessary?
n = nrow(train.set.multiple)
rmse.train = sqrt(1/n * sum(model.multiple$residuals^2))
n = nrow(test.set.multiple)
residuals = test.set.multiple$IC50.value.cisplatin - predict.multiple
rmse.test = sqrt(1/n * sum(residuals^2))
rmse.train
rmse.test

# model with specific biomarkers (GMDS, LRBA, ATXN1)
#same training and test set is used
model.multiple.biomarkers = lm(IC50.value.cisplatin ~ GMDS + LRBA + ATXN1, data = train.set.multiple)
summary(model.multiple.biomarkers)
plot(model.multiple.biomarkers)
predict.biomarkers = predict(model.multiple.biomarkers, newdata = test.set.multiple)

# for a better model, pca can be used
pca = prcomp(rges.ic50.biomarkers[, -2])
barplot(pca$rotation[, 1], horiz = TRUE, main = "PC1", col = "red")

# compute model with pcas instead of original variables
model.pca = lm(rges.ic50.biomarkers$IC50.value.cisplatin ~ pca$x)
summary(model.pca)
plot(model.pca)

