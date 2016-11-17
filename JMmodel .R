library("MASS")
library("nlme")
library("splines")
library("survival")
library("JM")
library("lattice")


##Joint Model 

## Data of patient 7
ND <- aids[aids$patient %in% c("7"), ] 

##Function 
pred.surv <- function(patient.data, prediction.time, fit.JM) {
  predJM <- survfitJM(fit.JM, newdata = patient.data, idVar = "patient", last.time = "Time", survTimes = prediction.time)
  return(predJM$summaries[[1]])
}


################################################ Model 1 ####################################################

##linear mixed effects model sqrt(cd4)~obstime & obstime:drug
fitLME <- lme(sqrt(CD4) ~ obstime + obstime:drug, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

#JM model 1* 
fit.JM11 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "weibull-AFT-GH")
fit.JM12 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "weibull-PH-GH")
fit.JM13 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "piecewise-PH-GH")
fit.JM14 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "Cox-PH-GH")
fit.JM15 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "spline-PH-GH")
fit.JM16 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM11)
summary(fit.JM12)
summary(fit.JM13)
summary(fit.JM14)
summary(fit.JM15)
summary(fit.JM16)




################################################ Model 2 ####################################################

##linear mixed effects model sqrt(cd4)~obstime & obstime:gender
fitLME2 <- lme(sqrt(CD4) ~ obstime + obstime:gender, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV2 <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)

##JM model 2*
fit.JM21 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "weibull-AFT-GH")
fit.JM22 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "weibull-PH-GH")
fit.JM23 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "piecewise-PH-GH")
fit.JM24 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "Cox-PH-GH")
fit.JM25 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "spline-PH-GH")
fit.JM26 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "ch-Laplace")

summary(fit.JM21)
summary(fit.JM22)
summary(fit.JM23)
summary(fit.JM24)
summary(fit.JM25)
summary(fit.JM26)



################################################ Model 3 ###################################################

##linear mixed effects model Time~obstime & CD4:gender
fitLME3 <- lme(sqrt(CD4) ~ obstime + CD4:gender, random = ~obstime | patient , data = aids)

##Cox model surv(time, death)~drug
fitSURV3 <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)

#JM model 3*
fit.JM31 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "weibull-AFT-GH")
fit.JM32 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "weibull-PH-GH")
fit.JM33 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "piecewise-PH-GH")
fit.JM34 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "Cox-PH-GH")
fit.JM35 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "spline-PH-GH")
fit.JM36 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "ch-Laplace")

summary(fit.JM31)
summary(fit.JM32)
summary(fit.JM33)
summary(fit.JM34)
summary(fit.JM35)
summary(fit.JM36)


################################################ Model 4 ###################################################

##linear mixed effects model Time~CD4 + prevOI + drug + AZT + start + stop + event
fitLME4 <- lme(sqrt(CD4) ~ CD4^2 + prevOI + drug + AZT +obstime^2, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV4 <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)

##JM model 4*
fit.JM41 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "weibull-AFT-GH")
fit.JM42 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "weibull-PH-GH")
fit.JM43 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "piecewise-PH-GH")
fit.JM44 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "Cox-PH-GH")
fit.JM45 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "spline-PH-GH")
fit.JM46 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "ch-Laplace")

summary(fit.JM41)
summary(fit.JM42)
summary(fit.JM43)
summary(fit.JM44)
summary(fit.JM45)
summary(fit.JM46)


################################################ Model 5 ###################################################

##linear mixed effects model sqrt(cd4)~obstime^3 & obstime:drug
fitLME5 <- lme(sqrt(CD4) ~ obstime^3 + obstime:drug, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV5 <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

##JM model 5*
fit.JM51 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "weibull-AFT-GH")
fit.JM52 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "weibull-PH-GH")
fit.JM53 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "piecewise-PH-GH")
fit.JM54 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "Cox-PH-GH")
fit.JM55 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "spline-PH-GH")
fit.JM56 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "ch-Laplace")

summary(fit.JM51)
summary(fit.JM52)
summary(fit.JM53)
summary(fit.JM54)
summary(fit.JM55)
summary(fit.JM56)

################################################ Model 6 ###################################################

##linear mixed effects model Time~obstime & CD4:gender
fitLME6 <- lme(sqrt(CD4) ~ obstime^2 + CD4:gender, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV6 <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)


##JM model 6*
fit.JM61 <- jointModel(fitLME6, fitSURV6, timeVar = "obstime", method = "weibull-AFT-GH")
fit.JM62 <- jointModel(fitLME6, fitSURV6, timeVar = "obstime", method = "weibull-PH-GH")
fit.JM63 <- jointModel(fitLME6, fitSURV6, timeVar = "obstime", method = "piecewise-PH-GH")
fit.JM64 <- jointModel(fitLME6, fitSURV6, timeVar = "obstime", method = "Cox-PH-GH")
fit.JM65 <- jointModel(fitLME6, fitSURV6, timeVar = "obstime", method = "spline-PH-GH")
fit.JM66 <- jointModel(fitLME6, fitSURV6, timeVar = "obstime", method = "ch-Laplace")

summary(fit.JM61)
summary(fit.JM62)
summary(fit.JM63)
summary(fit.JM64)
summary(fit.JM65)
summary(fit.JM66)



