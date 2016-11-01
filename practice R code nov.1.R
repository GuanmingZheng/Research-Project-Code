##Practice for r code
# load packages 'JM' and 'lattice'
library("MASS")
library("nlme")
library("splines")
library("survival")
library("JM")
library("lattice")

# Joint Model #


################################################## fit1 ####################################################### 

##linear mixed effects model sqrt(cd4)~obstime & obstime:drug
fitLME <- lme(sqrt(CD4) ~ obstime + obstime:drug, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)


############### fit1 the joint model with different method ##############

##method "weibull-AFT-GH"
fit.JM11 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM11)

##method "weibull-PH-GH"
fit.JM12 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM12)

##method "piecewise-PH-GH"
fit.JM13 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM13)

##method "Cox-PH-GH"
fit.JM14 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM14)

##method "spline-PH-GH" 
fit.JM15 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM15)

##method "ch-Laplace"
fit.JM16 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM16)







##################################################fit2####################################################### 
##linear mixed effects model sqrt(cd4)~obstime & obstime:gender
fitLME2 <- lme(sqrt(CD4) ~ obstime + obstime:gender, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV2 <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)

############### fit2 the joint model with different method ##############

##method "weibull-AFT-GH"
fit.JM21 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM21)

##method "weibull-PH-GH"
fit.JM22 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM22)

##method "piecewise-PH-GH"
fit.JM23 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM23)


##method "Cox-PH-GH"
fit.JM24 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM24)

##method "spline-PH-GH" 
fit.JM25 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM25)

##method "ch-Laplace"
fit.JM26 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM26)




##################################################fit3####################################################### 
##linear mixed effects model Time~obstime & CD4:gender
fitLME3 <- lme(sqrt(CD4) ~ obstime + CD4:gender, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV3 <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)

############### fit3 the joint model with different method ##############

##method "weibull-AFT-GH"
fit.JM31 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM31)

##method "weibull-PH-GH"
fit.JM32 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM32)

##method "piecewise-PH-GH"
fit.JM33 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM33)


##method "Cox-PH-GH"
fit.JM34 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM34)

##method "spline-PH-GH" 
fit.JM35 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM35)

##method "ch-Laplace"
fit.JM36 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM36)


##################################################fit4####################################################### 
##linear mixed effects model Time~CD4 & prevOI:gender
fitLME4 <- lme(sqrt(CD4) ~ Time + gender, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV4 <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)

############### fit4 the joint model with different method ##############

##method "weibull-AFT-GH"
fit.JM41 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM41)

##method "weibull-PH-GH"
fit.JM42 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM42)

##method "piecewise-PH-GH"
fit.JM43 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM43)


##method "Cox-PH-GH"
fit.JM44 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM44)

##method "spline-PH-GH" 
fit.JM45 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM45)

##method "ch-Laplace"
fit.JM46 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM46)



##################################################fit5####################################################### 
##linear mixed effects model Time~CD4 + prevOI + drug + AZT + start + stop + event
fitLME5 <- lme(sqrt(CD4) ~ CD4 + prevOI + drug + AZT +obstime, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV5 <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)

############### fit5 the joint model with different method ##############

##method "weibull-AFT-GH"
fit.JM51 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM51)

##method "weibull-PH-GH"
fit.JM52 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM52)

##method "piecewise-PH-GH"
fit.JM53 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM53)


##method "Cox-PH-GH"
fit.JM54 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM54)

##method "spline-PH-GH" 
fit.JM55 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM55)

##method "ch-Laplace"
fit.JM56 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM56)






##################################################fit7####################################################### with fit5
##linear mixed effects model Time~CD4 + prevOI + drug + AZT + start + stop + event
fitLME7 <- lme(sqrt(CD4) ~ CD4^2 + prevOI + drug + AZT +obstime^2, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV7 <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)
############### fit7 the joint model with different method ##############

##method "weibull-AFT-GH"
fit.JM71 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM71)

##method "weibull-PH-GH"
fit.JM72 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM72)

##method "piecewise-PH-GH"
fit.JM73 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM73)


##method "Cox-PH-GH"
fit.JM74 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM74)

##method "spline-PH-GH" 
fit.JM75 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM75)

##method "ch-Laplace"
fit.JM76 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM76)



################################################## fit8 ####################################################### with fit1

##linear mixed effects model sqrt(cd4)~obstime^3 & obstime:drug
fitLME8 <- lme(sqrt(CD4) ~ obstime^3 + obstime:drug, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV8 <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)


############### fit8 the joint model with different method ##############

##method "weibull-AFT-GH"
fit.JM81 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM81)

##method "weibull-PH-GH"
fit.JM82 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM82)

##method "piecewise-PH-GH"
fit.JM83 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM83)

##method "Cox-PH-GH"
fit.JM84 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM84)

##method "spline-PH-GH" 
fit.JM85 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM85)

##method "ch-Laplace"
fit.JM86 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM86)


################################################## fit9 ####################################################### with fit1

##linear mixed effects model sqrt(cd4)~obstime^3 & obstime:drug
fitLME9 <- lme(sqrt(CD4) ~ sqrt(obstime) + obstime:drug, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV9 <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)


############### fit9 the joint model with different method ##############

##method "weibull-AFT-GH"
fit.JM91 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM91)

##method "weibull-PH-GH"
fit.JM92 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM92)

##method "piecewise-PH-GH"
fit.JM93 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM93)

##method "Cox-PH-GH"
fit.JM94 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM94)

##method "spline-PH-GH" 
fit.JM95 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM95)

##method "ch-Laplace"
fit.JM96 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM96)

##################################################fit10####################################################### fit3
##linear mixed effects model Time~obstime & CD4:gender
fitLME10 <- lme(sqrt(CD4) ~ obstime^2 + CD4:gender, random = ~ obstime | patient, data = aids)

##Cox model surv(time, death)~drug
fitSURV10 <- coxph(Surv(Time, death) ~ gender, data = aids.id, x = TRUE)

############### fit3 the joint model with different method ##############

##method "weibull-AFT-GH"
fit.JM101 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM31)

##method "weibull-PH-GH"
fit.JM102 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM102)

##method "piecewise-PH-GH"
fit.JM103 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM103)


##method "Cox-PH-GH"
fit.JM104 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM104)

##method "spline-PH-GH" 
fit.JM105 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM105)

##method "ch-Laplace"
fit.JM106 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM106)




