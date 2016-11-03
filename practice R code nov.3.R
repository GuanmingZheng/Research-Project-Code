ND <- aids[aids$patient %in% c("7"), ]

ND.test = ND[,-c(7:12)]
Nd.test3 = ND[,-c(6:12)]
set.seed(123)
predSurv <- survfitJM(fit.JM11, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurv.test <- survfitJM(fit.JM11, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurv
predSurv.test
predSurv.test$summaries$'7'

plot(predSurv, which = "7", conf.int = TRUE)

### Give you patient information (survived up to time T)
### Give you a particular time (survival past time T+s)

data = ND.test
T.star = 16

set.seed(123)
predSurv.test <- survfitJM(fit.JM11, newdata = ND.test, idVar = "patient", last.time = "Time", survTimes = T.star)

pred.surv <- function(patient.data, prediction.time, fit.JM) {
  predJM <- survfitJM(fit.JM, newdata = patient.data, idVar = "patient", last.time = "Time", survTimes = prediction.time)
  return(predJM$summaries[[1]])
}

pred.surv(ND,16.5,fit.JM22)

##-----------------------------------------------Fit 1---------------------------------------------------##
##--------------fit1 JM11----------------## ##------"weibull-AFT-GH"------##

data= ND
ND.test = ND[,-c(7:12)]

set.seed(123)
predSurv <- survfitJM(fit.JM11, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurv.test <- survfitJM(fit.JM11, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurv
predSurv.test
plot(predSurv, which = "7", conf.int = TRUE)
##---------fit1 JM12----------## ##----------weibull-PH-GH-----------##

data= ND
ND.test = ND[,-c(7:12)]
set.seed(123)
predSurvJM12 <- survfitJM(fit.JM12, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM12.test <- survfitJM(fit.JM12, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM12
predSurvJM12.test
predSurvJM12.test$summaries$'7'
plot(predSurvJM12, which = "7", conf.int = TRUE)


##---------fit1 JM13----------## ##-----------piecewise-PH-GH---------##

data= ND
ND.test = ND[,-c(7:12)]
set.seed(123)
predSurvJM13 <- survfitJM(fit.JM13, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM13.test <- survfitJM(fit.JM13, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM13
predSurvJM13.test

plot(predSurvJM13, which = "7", conf.int = TRUE)



##---------fit1 JM14----------## ##-----------Cox-PH-GH---------## ## CANNOT RUN CODE##

fit.JM14 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM14)

data= ND
ND.test = ND[,-c(7:12)]
set.seed(123)
predSurvJM14 <- survfitJM(fit.JM14, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM14.test <- survfitJM(fit.JM14, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM14
predSurvJM14.test

plot(predSurvJM14, which = "7", conf.int = TRUE)


##--------------fit1 JM15-----------------## 

data= ND
ND.test = ND[,-c(7:12)]
set.seed(123)
predSurvJM15 <- survfitJM(fit.JM15, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM15.test <- survfitJM(fit.JM15, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM15
predSurvJM15.test

plot(predSurvJM15, which = "7", conf.int = TRUE)

##---------------fit1 JM16-------------------------##
data= ND
ND.test = ND[,-c(7:12)]
set.seed(123)
predSurvJM16 <- survfitJM(fit.JM16, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM16.test <- survfitJM(fit.JM16, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM16
predSurvJM16.test

plot(predSurvJM16, which = "7", conf.int = TRUE)


##----------------------------------------fit2---------------------------------------------##

##--------------fit2 JM21----------------## ##------"weibull-AFT-GH"------##

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM21 <- survfitJM(fit.JM21, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM21.test <- survfitJM(fit.JM21, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM21
predSurvJM21.test
plot(predSurv, which = "7", conf.int = TRUE)

##---------fit2 JM22----------## ##----------weibull-PH-GH-----------##
fit.JM22 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM22)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM22 <- survfitJM(fit.JM22, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM22.test <- survfitJM(fit.JM22, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM22
predSurvJM22.test
plot(predSurvJM22, which = "7", conf.int = TRUE)


##---------fit2 JM23----------## ##----------piecewise-PH-GH-----------##
fit.JM23 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM23)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM23 <- survfitJM(fit.JM23, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM23.test <- survfitJM(fit.JM23, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM23
predSurvJM23.test
plot(predSurvJM23, which = "7", conf.int = TRUE)

##---------fit2 JM24----------## ##----------Cox-PH-GH-----------##
fit.JM24 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM24)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM24 <- survfitJM(fit.JM24, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM24.test <- survfitJM(fit.JM24, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM24
predSurvJM24.test
plot(predSurvJM24, which = "7", conf.int = TRUE)



##---------fit2 JM25----------## ##----------spline-PH-GH-----------##
fit.JM25 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM25)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM25 <- survfitJM(fit.JM25, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM25.test <- survfitJM(fit.JM25, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM25
predSurvJM25.test
plot(predSurvJM25, which = "7", conf.int = TRUE)


##---------fit2 JM26----------## ##----------ch-Laplace-----------##
fit.JM26 <- jointModel(fitLME2, fitSURV2, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM26)
data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM26 <- survfitJM(fit.JM26, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM26.test <- survfitJM(fit.JM26, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM26
predSurvJM26.test
plot(predSurvJM26, which = "7", conf.int = TRUE)




##----------------------------------------fit3---------------------------------------------##

##--------------fit3 JM31----------------## ##------"weibull-AFT-GH"------##

fit.JM31 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM31)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM31 <- survfitJM(fit.JM31, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM31.test <- survfitJM(fit.JM31, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM31
predSurvJM31.test
plot(predSurvJM31, which = "7", conf.int = TRUE)


##--------------fit3 JM32----------------## ##------"weibull-PH-GH"------##

fit.JM32 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM32)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM32 <- survfitJM(fit.JM32, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM32.test <- survfitJM(fit.JM32, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM32
predSurvJM32.test
plot(predSurvJM32, which = "7", conf.int = TRUE)

##--------------fit3 JM33----------------## ##------"weibull-PH-GH"------##

fit.JM32 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM32)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM33 <- survfitJM(fit.JM33, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM33.test <- survfitJM(fit.JM33, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM33
predSurvJM33.test
plot(predSurvJM33, which = "7", conf.int = TRUE)

##--------------fit3 JM34----------------## ##------"Cox-PH-GH"------##


fit.JM34 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM34)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM34 <- survfitJM(fit.JM34, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM34.test <- survfitJM(fit.JM34, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM34
predSurvJM34.test
plot(predSurvJM34, which = "7", conf.int = TRUE)


##--------------fit3 JM35----------------## ##------"spline-PH-GH"------##

fit.JM35 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM35)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM35 <- survfitJM(fit.JM35, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM35.test <- survfitJM(fit.JM35, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM35
predSurvJM35.test
plot(predSurvJM35, which = "7", conf.int = TRUE)


##--------------fit3 JM36----------------## ##------"ch-Laplace"------##

fit.JM36 <- jointModel(fitLME3, fitSURV3, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM36)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM36 <- survfitJM(fit.JM36, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM36.test <- survfitJM(fit.JM36, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM36
predSurvJM36.test
plot(predSurvJM36, which = "7", conf.int = TRUE)


##----------------------------------------fit4---------------------------------------------##

##--------------fit4 JM41----------------## ##------"weibull-AFT-GH"------##


fit.JM41 <- jointModel(fitLME4, fitSURV4, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM41)

data= ND
ND.test = ND[,-c(8:12)]

set.seed(123)
predSurvJM41 <- survfitJM(fit.JM41, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM41.test <- survfitJM(fit.JM41, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM41
predSurvJM41.test
plot(predSurvJM41, which = "7", conf.int = TRUE)






##----------------------------------------fit5---------------------------------------------##

##--------------fit5 JM51----------------## ##------"weibull-AFT-GH"------##

fit.JM51 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM51)

data5= ND
ND.test5 = ND[,-c(10:12)]

set.seed(123)
predSurvJM51 <- survfitJM(fit.JM51, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM51.test <- survfitJM(fit.JM51, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM51
predSurvJM51.test
plot(predSurvJM51, which = "7", conf.int = TRUE)


##--------------fit5 JM52----------------## ##------"weibull-PH-GH"------##

fit.JM52 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM52)

data5= ND
ND.test5 = ND[,-c(10:12)]

set.seed(123)
predSurvJM52 <- survfitJM(fit.JM52, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM52.test <- survfitJM(fit.JM52, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM52
predSurvJM52.test
plot(predSurvJM52, which = "7", conf.int = TRUE)

##--------------fit3 JM53----------------## ##------"weibull-PH-GH"------##

fit.JM53 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM53)

data5= ND
ND.test5 = ND[,-c(10:12)]

set.seed(123)
predSurvJM53 <- survfitJM(fit.JM53, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM53.test <- survfitJM(fit.JM53, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM53
predSurvJM53.test
plot(predSurvJM53, which = "7", conf.int = TRUE)

##--------------fit5 JM54----------------## ##------"Cox-PH-GH"------##

fit.JM54 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM54)

data5= ND
ND.test5 = ND[,-c(10:12)]

set.seed(123)
predSurvJM54 <- survfitJM(fit.JM54, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM54.test <- survfitJM(fit.JM54, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM54
predSurvJM54.test
plot(predSurvJM54, which = "7", conf.int = TRUE)


##--------------fit5 JM55----------------## ##------"spline-PH-GH"------##

fit.JM55 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM55)

data5= ND
ND.test5 = ND[,-c(10:12)]

set.seed(123)
predSurvJM55 <- survfitJM(fit.JM55, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM55.test <- survfitJM(fit.JM55, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM55
predSurvJM55.test
plot(predSurvJM55, which = "7", conf.int = TRUE)


##--------------fit5 JM56----------------## ##------"ch-Laplace"------##

fit.JM56 <- jointModel(fitLME5, fitSURV5, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM56)

data5= ND
ND.test5 = ND[,-c(10:12)]

set.seed(123)
predSurvJM56 <- survfitJM(fit.JM56, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM56.test <- survfitJM(fit.JM56, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM56
predSurvJM56.test
plot(predSurvJM56, which = "7", conf.int = TRUE)



##----------------------------------------fit7---------------------------------------------##

##--------------fit7 JM71----------------## ##------"weibull-AFT-GH"------##

fit.JM71 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM71)

data7= ND
ND.test7 = ND[,-c(10:12)]

set.seed(123)
predSurvJM71 <- survfitJM(fit.JM71, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM71.test <- survfitJM(fit.JM71, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM71
predSurvJM71.test
plot(predSurvJM71, which = "7", conf.int = TRUE)


##--------------fit5 JM72----------------## ##------"weibull-PH-GH"------##

fit.JM72 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM72)

data7= ND
ND.test7 = ND[,-c(10:12)]

set.seed(123)
predSurvJM72 <- survfitJM(fit.JM72, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM72.test <- survfitJM(fit.JM72, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM72
predSurvJM72.test
plot(predSurvJM72, which = "7", conf.int = TRUE)

##--------------fit7 JM73----------------## ##------"weibull-PH-GH"------##

fit.JM73 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM73)


data7= ND
ND.test7 = ND[,-c(10:12)]

set.seed(123)
predSurvJM73 <- survfitJM(fit.JM73, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM73.test <- survfitJM(fit.JM73, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM73
predSurvJM73.test
plot(predSurvJM73, which = "7", conf.int = TRUE)

##--------------fit7 JM74----------------## ##------"Cox-PH-GH"------##

fit.JM74 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM74)

data7= ND
ND.test7 = ND[,-c(10:12)]

set.seed(123)
predSurvJM74 <- survfitJM(fit.JM74, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM74.test <- survfitJM(fit.JM74, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM74
predSurvJM74.test
plot(predSurvJM74, which = "7", conf.int = TRUE)


##--------------fit7 JM75----------------## ##------"spline-PH-GH"------##

fit.JM75 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM75)


data7= ND
ND.test7 = ND[,-c(10:12)]

set.seed(123)
predSurvJM75 <- survfitJM(fit.JM75, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM75.test <- survfitJM(fit.JM75, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM75
predSurvJM75.test
plot(predSurvJM75, which = "7", conf.int = TRUE)


##--------------fit7 JM76----------------## ##------"ch-Laplace"------##

fit.JM76 <- jointModel(fitLME7, fitSURV7, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM76)

data7= ND
ND.test7 = ND[,-c(10:12)]

set.seed(123)
predSurvJM76 <- survfitJM(fit.JM76, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM76.test <- survfitJM(fit.JM76, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM76
predSurvJM76.test
plot(predSurvJM76, which = "7", conf.int = TRUE)




##-----------------------------------------------Fit 8---------------------------------------------------##
##--------------fit8 JM81----------------## ##------"weibull-AFT-GH"------##

fit.JM81 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM81)

data8= ND
ND.test8 = ND[,-c(7:12)]

set.seed(123)
predSurvJM81 <- survfitJM(fit.JM81, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM81.test <- survfitJM(fit.JM81, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM81
predSurvJM81.test
plot(predSurvJM81, which = "7", conf.int = TRUE)

##---------fit8 JM82----------## ##----------weibull-PH-GH-----------##

fit.JM82 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM82)

data8= ND
ND.test8 = ND[,-c(7:12)]

set.seed(123)
predSurvJM82 <- survfitJM(fit.JM82, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM82.test <- survfitJM(fit.JM82, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM82
predSurvJM82.test
plot(predSurvJM82, which = "7", conf.int = TRUE)



##---------fit8 JM83----------## ##-----------piecewise-PH-GH---------##

fit.JM83 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM83)

data8= ND
ND.test8 = ND[,-c(7:12)]

set.seed(123)
predSurvJM83 <- survfitJM(fit.JM83, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM83.test <- survfitJM(fit.JM83, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM83
predSurvJM83.test

plot(predSurvJM83, which = "7", conf.int = TRUE)



##---------fit8 JM84----------## ##-----------Cox-PH-GH---------## ## CANNOT RUN CODE##

fit.JM84 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM84)

data8= ND
ND.test8 = ND[,-c(7:12)]

set.seed(123)
predSurvJM84 <- survfitJM(fit.JM84, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM84.test <- survfitJM(fit.JM84, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM84
predSurvJM84.test

plot(predSurvJM84, which = "7", conf.int = TRUE)


##--------------fit8 JM85-----------------## 

fit.JM85 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM85)

data8= ND
ND.test8 = ND[,-c(7:12)]

set.seed(123)
predSurvJM85 <- survfitJM(fit.JM85, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM85.test <- survfitJM(fit.JM85, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM85
predSurvJM85.test

plot(predSurvJM85, which = "7", conf.int = TRUE)

##---------------fit8 JM86---------## ----------"ch-Laplace"------##

fit.JM86 <- jointModel(fitLME8, fitSURV8, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM86)

data8= ND
ND.test8 = ND[,-c(7:12)]
set.seed(123)
predSurvJM86 <- survfitJM(fit.JM86, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM86.test <- survfitJM(fit.JM86, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM86
predSurvJM86.test

plot(predSurvJM86, which = "7", conf.int = TRUE)



##-----------------------------------------------Fit 9---------------------------------------------------##
##--------------fit JM91----------------## ##------"weibull-AFT-GH"------##

fit.JM91 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM91)

data9= ND
ND.test9 = ND[,-c(7:12)]

set.seed(123)
predSurvJM91 <- survfitJM(fit.JM91, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM91.test <- survfitJM(fit.JM91, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM91
predSurvJM91.test
plot(predSurvJM91, which = "7", conf.int = TRUE)

##---------fit9 JM92----------## ##----------weibull-PH-GH-----------##

fit.JM92 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM92)

data9= ND
ND.test9 = ND[,-c(7:12)]

set.seed(123)
predSurvJM92 <- survfitJM(fit.JM92, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM92.test <- survfitJM(fit.JM92, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM92
predSurvJM92.test
plot(predSurvJM92, which = "7", conf.int = TRUE)



##---------fit9 JM93----------## ##-----------piecewise-PH-GH---------##

fit.JM93 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM93)

data9= ND
ND.test9 = ND[,-c(7:12)]

set.seed(123)
predSurvJM93 <- survfitJM(fit.JM93, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM93.test <- survfitJM(fit.JM93, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM93
predSurvJM93.test

plot(predSurvJM93, which = "7", conf.int = TRUE)



##---------fit9 JM94----------## ##-----------Cox-PH-GH---------## ## CANNOT RUN CODE##

fit.JM94 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM94)

data9= ND
ND.test9 = ND[,-c(7:12)]

set.seed(123)
predSurvJM94 <- survfitJM(fit.JM94, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM94.test <- survfitJM(fit.JM94, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM94
predSurvJM94.test

plot(predSurvJM94, which = "7", conf.int = TRUE)


##--------------fit9 JM95-----------------## 

fit.JM95 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM95)

data9= ND
ND.test9 = ND[,-c(7:12)]

set.seed(123)
predSurvJM95 <- survfitJM(fit.JM95, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM95.test <- survfitJM(fit.JM95, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM95
predSurvJM95.test

plot(predSurvJM95, which = "7", conf.int = TRUE)

##---------------fit9 JM96---------## ----------"ch-Laplace"------##

fit.JM96 <- jointModel(fitLME9, fitSURV9, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM96)

data9= ND
ND.test9 = ND[,-c(7:12)]
set.seed(123)
predSurvJM96 <- survfitJM(fit.JM96, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM96.test <- survfitJM(fit.JM96, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM96
predSurvJM96.test

plot(predSurvJM96, which = "7", conf.int = TRUE)


##----------------------------------------fit10---------------------------------------------##

##--------------fit10 JM101----------------## ##------"weibull-AFT-GH"------##

fit.JM101 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "weibull-AFT-GH")
summary(fit.JM31)

data10= ND
ND.test10 = ND[,-c(8:12)]

set.seed(123)
predSurvJM101 <- survfitJM(fit.JM101, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM101.test <- survfitJM(fit.JM101, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM101
predSurvJM101.test
plot(predSurvJM101, which = "7", conf.int = TRUE)


##--------------fit10 JM102----------------## ##------"weibull-PH-GH"------##

fit.JM102 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "weibull-PH-GH")
summary(fit.JM102)
data10= ND
ND.test10 = ND[,-c(8:12)]

set.seed(123)
predSurvJM102 <- survfitJM(fit.JM102, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM102.test <- survfitJM(fit.JM102, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM102
predSurvJM102.test
plot(predSurvJM102, which = "7", conf.int = TRUE)

##--------------fit10 JM103----------------## ##------"weibull-PH-GH"------##

fit.JM103 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM103)

data10= ND
ND.test10 = ND[,-c(8:12)]

set.seed(123)
predSurvJM103 <- survfitJM(fit.JM103, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM103.test <- survfitJM(fit.JM103, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM103
predSurvJM103.test
plot(predSurvJM103, which = "7", conf.int = TRUE)

##--------------fit10 JM104----------------## ##------"Cox-PH-GH"------##


fit.JM104 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "Cox-PH-GH")
summary(fit.JM104)

data10= ND
ND.test10 = ND[,-c(8:12)]

set.seed(123)
predSurvJM104 <- survfitJM(fit.JM104, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM104.test <- survfitJM(fit.JM104, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM104
predSurvJM104.test
plot(predSurvJM104, which = "7", conf.int = TRUE)


##--------------fit10 JM105----------------## ##------"spline-PH-GH"------##

fit.JM105 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "spline-PH-GH")
summary(fit.JM105)

data10= ND
ND.test10 = ND[,-c(8:12)]

set.seed(123)
predSurvJM105 <- survfitJM(fit.JM105, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM105.test <- survfitJM(fit.JM105, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM105
predSurvJM105.test
plot(predSurvJM105, which = "7", conf.int = TRUE)


##--------------fit3 JM36----------------## ##------"ch-Laplace"------##

fit.JM106 <- jointModel(fitLME10, fitSURV10, timeVar = "obstime", method = "ch-Laplace")
summary(fit.JM106)

data10= ND
ND.test10 = ND[,-c(8:12)]

set.seed(123)
predSurvJM106 <- survfitJM(fit.JM106, newdata = ND, idVar = "patient", last.time = "Time")

set.seed(123)
predSurvJM106.test <- survfitJM(fit.JM106, newdata = ND.test, idVar = "patient", last.time = "Time")

predSurvJM106
predSurvJM106.test
plot(predSurvJM106, which = "7", conf.int = TRUE)




