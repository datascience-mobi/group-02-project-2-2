# limit the data set to numerical values
rges.ic50 = drug_activity_rges[,c(2,9)]

# compute correlation
cor(rges.ic50$RGES, rges.ic50$IC50.value)
# check significance
cor.test(rges.ic50$RGES, rges.ic50$IC50.value)

## take 200 random samples to form the training set
i.train = sample(1:nrow(rges.ic50), 200)
## 
rges.ic50.train = rges.ic50[i.train, ]
rges.ic50.test = rges.ic50[-i.train, ]
#learn lm on training set
l.train = lm(IC50.value ~ RGES, data = rges.ic50.train)
summary(l.train)

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


