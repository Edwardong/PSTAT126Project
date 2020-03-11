library(base)
library(car)
library(carData)
library(corrplot)
library(datasets)
library(forcats)
library(ggplot2)
library(graphics)
library(methods)
library(readr)
library(dplyr)
library(leaps)
library(purrr)
library(MASS)
library(tibble)
library(tidyr)
library(utils)

#load file into datagram
QSAR_fish <- read.csv("qsar_fish_toxicity.csv")
names(QSAR_fish) <- c("CIC0", "SM1_Dz", "GATS1i", "NdsCH", "NdssC", "MLOGP","LC50")

#initialize variables
CIC0 <- QSAR_fish$CIC0
SM1_Dz <- QSAR_fish$SM1_Dz
GATS1i <- QSAR_fish$GATS1i
NdsCH <- QSAR_fish$NdsCH
NdssC <- QSAR_fish$NdssC
MLOGP <- QSAR_fish$MLOGP

#Response
LC50 <- QSAR_fish$LC50

#Scatterplot
pairs(LC50~CIC0+SM1_Dz+GATS1i+NdsCH+NdssC+MLOGP,QSAR_fish,main = "Scatterplot Matrix of QSAR_fish data")
#pairs(LC50^(0.8585859)~MLOGP)

#correlation
cor(subset(QSAR_fish,select = -c(LC50)))
#Variance Inflation Factor
#vif(subset(QSAR_fish,select = -c(LC50)))

#First model
full.mod <- lm(LC50~CIC0+SM1_Dz+GATS1i+factor(NdsCH)+factor(NdssC)+MLOGP)
#plots, observe Non-constant error variance
plot(full.mod)
#histgram
hist(resid(full.mod),main = "Histogram of the Residuals", xlab = "Residuals")

#plot(CIC0, resid(full.mod), main = "Residuals vs. CIC0", xlab = "CIC0", ylab = "residuals")

#plot(SM1_Dz, resid(full.mod), main = "Residuals vs. SM1_Dz", xlab = "SM1_Dz", ylab = "residuals")
#plot(GATS1i, resid(full.mod), main = "Residuals vs. GATS1i", xlab = "GATS1i", ylab = "residuals")
#plot(factor(NdsCH), resid(full.mod), main = "Residuals vs. NdsCH", xlab = "NdsCH", ylab = "residuals")
#plot(factor(NdssC), resid(full.mod), main = "Residuals vs. NdssC", xlab = "NdssC", ylab = "residuals")
#plot(MLOGP, resid(full.mod), main = "Residuals vs. MLOGP", xlab = "MLOGP", ylab = "residuals")

#boccox
potential.trans <- boxcox(full.mod, lambda=seq(-1,1,length=10))
potential.trans$x[which.max(potential.trans$y)]

#Second model dealt with Non-constant error variance with Transformation on Response
full.mod2 <-lm(LC50^(0.8585859)~CIC0+SM1_Dz+GATS1i+factor(NdsCH)+factor(NdssC)+MLOGP)
plot(full.mod2)
hist(resid(full.mod2),main = "Histogram of the Residuals", xlab = "Residuals")

#Stepwise Regression using AIC
mod0 = lm(LC50^(0.8585859)~1)
mod.upper = lm(LC50^(0.8585859)~CIC0+SM1_Dz+GATS1i+NdsCH+NdssC+MLOGP)
summary(mod.upper)
step(mod0,scope = list(lower=mod0,upper=mod.upper))

#Best subset regression
mod.bsubsets <- regsubsets(cbind(CIC0,SM1_Dz,GATS1i,NdsCH,NdssC,MLOGP),LC50^(0.8585859))
summary.bsubsets <- summary(mod.bsubsets)
#with adjusted R^2 criteria
summary.bsubsets$adjr2
summary.bsubsets$which
#Mallows CP
summary.bsubsets$cp

#Stepwise regression using F-tests
add1(mod0,~.+CIC0+SM1_Dz+GATS1i+NdsCH+NdssC+MLOGP,test="F")
mod1 = update(mod0, ~.+MLOGP)
add1(mod1,~.+CIC0+SM1_Dz+GATS1i+NdsCH+NdssC,test="F")
mod2 = update(mod1, ~.+SM1_Dz)
summary(mod2)
add1(mod2,~.+CIC0+GATS1i+NdsCH+NdssC,test="F")
mod3 = update(mod2, ~.+NdsCH)
summary(mod3)
add1(mod3,~.+CIC0+GATS1i+NdssC,test="F")
mod4 = update(mod3, ~.+GATS1i)
summary(mod4)
add1(mod4,~.+CIC0+NdssC,test="F")
mod5 = update(mod4, ~.+CIC0)
summary(mod5)
add1(mod5,~.+NdssC,test="F")
#not adding NdssC

#Delete influential points using Studentized Deleted Residuals
which(abs(rstudent(full.mod2)) > 3)
Data_deleted <- QSAR_fish[-c(66,121,170,183,220,238,260,288,373,448,681),]
CIC0_d <- Data_deleted$CIC0
SM1_Dz_d <- Data_deleted$SM1_Dz
GATS1i_d <- Data_deleted$GATS1i
NdsCH_d <- Data_deleted$NdsCH
NdssC_d <- Data_deleted$NdssC
MLOGP_d <- Data_deleted$MLOGP

LC50_d <- Data_deleted$LC50

full.mod4 <- lm(LC50_d^(0.8585859)~MLOGP_d+SM1_Dz_d+factor(NdsCH_d)+GATS1i_d+CIC0_d)

#deleted base on DFFITS value
dffits(full.mod2)
dim(QSAR_fish)
num = 2 * sqrt((7 + 1)/(908-7-1))
which(dffits(full.mod2) > num)
Data_deleted_2 <- QSAR_fish[-c(19,31,45,62,66,73,75,76,77,78,84,85,121,152,166,170,178,183,220,231,238,250,255,257,260,282,288,303,348,350,364,373,417,418,429,443,448,455,463,483,503,505,551,552,563,578,681,692,695,709,724,726,761,768,784,826,890,893,901,905,906),]
CIC0_d2 <- Data_deleted_2$CIC0
SM1_Dz_d2 <- Data_deleted_2$SM1_Dz
GATS1i_d2 <- Data_deleted_2$GATS1i
NdsCH_d2 <- Data_deleted_2$NdsCH
NdssC_d2 <- Data_deleted_2$NdssC
MLOGP_d2 <- Data_deleted_2$MLOGP

LC50_d2 <- Data_deleted_2$LC50
full.mod5 <- lm(LC50_d2^(0.8585859)~MLOGP_d2+SM1_Dz_d2+factor(NdsCH_d2)+GATS1i_d2+CIC0_d2)

#Question1
reduced.mod <- lm(LC50_d2^(0.8585859)~1)
full.mod6 <- lm(LC50_d2^(0.8585859)~MLOGP_d2+SM1_Dz_d2+factor(NdsCH_d2)+GATS1i_d2+CIC0_d2)
anova(reduced.mod, full.mod6)

#Question2
new.df = data.frame(MLOGP_d2 = mean(MLOGP_d2), SM1_Dz_d2 = mean(SM1_Dz_d2), CIC0_d2 = 3, GATS1i_d2=mean(GATS1i_d2), NdsCH_d2 = 0,NdssC_d2 = 0)
predicted = predict(full.mod8, new.df, se.fit = TRUE, interval = "confidence", level = 0.95, type = "response")
predicted$fit

#Question3
interaction_mod1 = update(full.mod5, ~. + CIC0_d3*MLOGP_d3)
plot(CIC0_d3*MLOGP_d3, resid(full.mod5),main = "Residuals vs. CIC0 * MLOGP", xlab = "CIC0 * MLOGP",
     ylab = "Residuals")
anova(full.mod5,interaction_mod1)

interaction_mod2 = update(interaction_mod1, ~. + CIC0_d3*SM1_Dz_d3)
plot(CIC0_d3*SM1_Dz_d3, resid(interaction_mod1),main = "Residuals vs. CIC0 * SM1_Dz", xlab = "CIC0 * SM1_Dz",
     ylab = "Residuals")
anova(interaction_mod1,interaction_mod2)

interaction_mod2.1 = update(interaction_mod1, ~. + CIC0_d3*GATS1i_d3)
plot(CIC0_d3*GATS1i_d3, resid(interaction_mod1),main = "Residuals vs. CIC0 * GATS1i_d3", xlab = "CIC0 * GATS1i_d3",
     ylab = "Residuals")
anova(interaction_mod1,interaction_mod2.1)

interaction_mod3 = update(interaction_mod2, ~. + CIC0_d3*GATS1i_d3)
plot(CIC0_d3*GATS1i_d3, resid(interaction_mod2),main = "Residuals vs. Adjusted model with CIC0 * GATS1i", xlab = "CIC0 * GATS1i",
     ylab = "Residuals")
anova(interaction_mod2,interaction_mod3)

plot(interaction_mod2)

hist(resid(interaction_mod2), main = "Histogram of the Residuals", xlab = "Residuals")

summary(interaction_mod2)


