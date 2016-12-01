
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

## Functions
set.Ns<-function(data,s){
  potential.people <-data$patient[data$Time >= s]
  return(unique(potential.people))
}

Surv.pred<-function(patient.data, s, t, fit.JM){
  data1 <- patient.data
  data1$Time <- s
  sub.data1 <- subset(data1,data1$obstime < s+t)
  predJM <- survfitJM(fit.JM, newdata = sub.data1, idVar = "patient", last.time = "Time", survTimes = s+t)
  return(predJM$summaries[[1]][2])
}

library(survival)
KM.cens <- survfit(Surv(Time, !death) ~ 1,  type="kaplan-meier", conf.type="log", data=aids.id)


hat.G <- function(KM.cens, t) {
  abs.diff = KM.cens$time-t
  return(KM.cens$surv[max(which(abs.diff<=0))])
}

hat.BS.st <- function(data, fit.JM, s, t) {
  BS.total = 0
  people.s = set.Ns(aids,s)
  
  for (patient in people.s) {
    patient.data = aids[aids$patient == patient,]
    
    patient.Shat = Surv.pred(patient.data, s, t, fit.JM)
    
    patient.y = as.numeric(patient.data$Time[1] > t + s)
    
    if(patient.y == 0) {
      censored = (patient.data$death[1] == 0)
      denominator = (hat.G(KM.cens, patient.data$Time[1])) / (hat.G(KM.cens, s)) 
      BS.total = BS.total + censored*(patient.y - patient.Shat)^2 / denominator
    } else {
      denominator = (hat.G(KM.cens, t+s)) / (hat.G(KM.cens, s)) 
      BS.total = BS.total + (patient.y - patient.Shat)^2 / denominator
    }
    print(patient)
    
  }
  
  return(BS.total/length(people.s))
}

