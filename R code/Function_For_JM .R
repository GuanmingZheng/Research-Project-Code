library("MASS")
library("nlme")
library("splines")
library("survival")
library("JM")
library("lattice")

##Data 
ND <- aids[aids$patient %in% c("7"), ]
ND.test = ND[,-c(7:12)]

##Model 
fitLME <- lme(sqrt(CD4) ~ obstime + obstime:drug, random = ~ obstime | patient, data = aids)
fitSURV <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)
fit.JM11 <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "weibull-AFT-GH")


##function 1 (pred.surv function)

pred.surv <- function(patient.data, prediction.time, fit.JM) {
  predJM <- survfitJM(fit.JM, newdata = patient.data, idVar = "patient", last.time = "Time", survTimes = prediction.time)
  return(predJM$summaries[[1]])
}


#Test
pred.surv(ND,16.5,fit.JM11)
pred.surv(ND,17,fit.JM11)
pred.surv(ND,17,fit.JM12)
pred.surv(ND,17,fit.JM32)




##fuction 2 (Quadratic error of prediction)

#Part1:function for Ns, which is the number of subjects still at risk at time s.

Ns<-function(patient.data,s,t){
  n<-patient.data$Time >= s 
  N<-length(which(n == TRUE))
  return(N) 
}


##Test##

Ns(aids,15,1)
Ns(aids,13,1)
Ns(aids,3,2)



#part2:
S.Pred<-function(patient.data,t,fit.JM){
  
  data1<-patient.data
  data1$Time = t
  data2<-subset(data1,data1$obstime >t)
  predJM <- survfitJM(fit.JM, newdata = data2, idVar = "patient", last.time = "Time", survTimes = t)
  return(predJM$summaries[[1]])
}




