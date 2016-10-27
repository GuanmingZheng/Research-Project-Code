# load packages 'JM' and 'lattice'
library("MASS")
library("nlme")
library("splines")
library("survival")
library("JM")
library("lattice")



# Descriptive Plots #


# longitudinal outcome
##GRPAH
##1##graph CD4~obstime and drug 

xyplot(sqrt(CD4) ~ obstime | drug, group = patient, data = aids, 
       xlab = "Months", ylab = expression(sqrt("CD4")), col = 1, type = "l")

##2##graph cd4~obstime and gender

xyplot(sqrt(CD4)~ obstime | gender, group = patient, data= aids, xlab = "times",ylab=expression(sqrt("CD4")),col=1,type="l")

##3##graph cd4~obstime and prevOI
xyplot(sqrt(CD4)~ obstime | prevOI, group = patient, data= aids, xlab = "times",ylab=expression(sqrt("CD4")),col=1,type="l")

##4##graph cd4~obstime and AZT
xyplot(sqrt(CD4)~ obstime | AZT, group = patient, data= aids, xlab = "times",ylab=expression(sqrt("CD4")),col=1,type="l")




# survival outcome
plot(survfit(Surv(Time, death) ~ drug, data = aids.id), conf.int = FALSE, 
     mark.time = TRUE, col = c("black", "red"), lty = 1:2, 
     ylab = "Survival", xlab = "Months") 
legend("topright", c("ddC", "ddI"), lty = 1:2, col = c("black", "red"), 
       bty = "n")


# Naive Cox Model #

##Cox model 1
td.Cox <- coxph(Surv(start, stop, event) ~ drug + sqrt(CD4), data = aids)
summary(td.Cox)

##Cox model 2
td.Cox2 <- coxph(Surv(start, stop, event) ~ gender + sqrt(CD4)^2, data = aids)
summary(td.Cox2)

##Cox model 3
td.Cox3 <- coxph(Surv(start, stop, event) ~ gender + drug*AZT + sqrt(CD4)^2, data = aids)
summary(td.Cox3)

##Cox model 4
td.Cox4 <- coxph(Surv(start, stop, event) ~ gender + drug + sqrt(CD4)^2, data = aids)
summary(td.Cox4)

##Cox model 5
td.Cox5 <- coxph(Surv(start, stop, event) ~ drug + gender^2 + sqrt(CD4)^2, data = aids)
summary(td.Cox5)

# Joint Model #


# first we fit a linear mixed effects model with lme(),
fitLME <- lme(sqrt(CD4) ~ obstime + obstime:drug, random = ~ obstime | patient, data = aids)
# and a Cox model with coxph() -- you need 'x = TRUE' in the call to coxph()
fitSURV <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# fit the joint model using function jointModel()
fit.JM <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM)


##fit2
fitLMEb <- lme(sqrt(CD4) ~ obstime + obstime:gender, random = ~ obstime | patient, data = aids)
fitSURVb <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)
fit.JMb <- jointModel(fitLMEb, fitSURVb, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JMb)

##fit3
fitLMEc <- lme(sqrt(CD4)^2 ~ Time + obstime:prevOI, random = ~ obstime | patient, data = aids)
fitSURVc <- coxph(Surv(Time, death) ~ prevOI, data = aids.id, x = TRUE)
fit.JMc <- jointModel(fitLMEc, fitSURVc, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JMc)

##fit4
fitLMEd <- lme(sqrt(CD4) ~ obstime + obstime:drug, random = ~ obstime | patient, data = aids)
fitSURVd <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)
fit.JMd <- jointModel(fitLMEd, fitSURVd, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JMd)




# likelihood ratio test for no treatment effect in the survival process:
# we fit the joint model under the null hypothesis
fitSURV2 <- coxph(Surv(Time, death) ~ 1, data = aids.id, x = TRUE)
fit.JM2 <- jointModel(fitLME, fitSURV2, timeVar = "obstime", method = "piecewise-PH-GH")
# the likelihood ratio test is performed with the anova() method
anova(fit.JM2, fit.JM)
