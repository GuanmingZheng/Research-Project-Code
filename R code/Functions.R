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
  people.s = set.Ns(aids,s) ##get the number of subject
  
  for (patient in people.s) {
    patient.data = aids[aids$patient == patient,] ## get the data set for patient 
    
    patient.Shat = Surv.pred(patient.data, s, t, fit.JM)
    
    patient.y = as.numeric(patient.data$Time[1] > t + s)  ## get a value to patient.y if the time of patient i is greater than t+s
    
    if(patient.y == 0) {
      censored = (patient.data$death[1] == 0) ##death = 0 means censored 
      denominator = (hat.G(KM.cens, patient.data$Time[1])) / (hat.G(KM.cens, s)) 
      BS.total = BS.total + censored*(patient.y - patient.Shat)^2 / denominator
    } else { 
      denominator = (hat.G(KM.cens, t+s)) / (hat.G(KM.cens, s)) 
      BS.total = BS.total + (patient.y - patient.Shat)^2 / denominator
    }
   ## print(patient)
 
 }
 
  return(BS.total/length(people.s))
}




############Draw plot  

S=seq(0,10,0.5)
S
t1=3
t2=5
z1=hat.BS.st(ND,fit.JM11,0,3)
z2=hat.BS.st(ND,fit.JM11,0.5,3)
z3=hat.BS.st(ND,fit.JM11,1,3)
z4=hat.BS.st(ND,fit.JM11,1.5,3)
z5=hat.BS.st(ND,fit.JM11,2,3)
z6=hat.BS.st(ND,fit.JM11,2.5,3)
z7=hat.BS.st(ND,fit.JM11,3,3)
z8=hat.BS.st(ND,fit.JM11,3.5,3)
z9=hat.BS.st(ND,fit.JM11,4,3)
z10=hat.BS.st(ND,fit.JM11,4.5,3)
z11=hat.BS.st(ND,fit.JM11,5,3)
z12=hat.BS.st(ND,fit.JM11,5.5,3)
z13=hat.BS.st(ND,fit.JM11,6,3)
z14=hat.BS.st(ND,fit.JM11,6.5,3)
z15=hat.BS.st(ND,fit.JM11,7,3)
z16=hat.BS.st(ND,fit.JM11,7.5,3)
z17=hat.BS.st(ND,fit.JM11,8,3)
z18=hat.BS.st(ND,fit.JM11,8.5,3)
z19=hat.BS.st(ND,fit.JM11,9,3)
z20=hat.BS.st(ND,fit.JM11,9.5,3)
z21=hat.BS.st(ND,fit.JM11,10,3)

q1=hat.BS.st(ND,fit.JM11,0,5)
q2=hat.BS.st(ND,fit.JM11,0.5,5)
q3=hat.BS.st(ND,fit.JM11,1,5)
q4=hat.BS.st(ND,fit.JM11,1.5,5)
q5=hat.BS.st(ND,fit.JM11,2,5)
q6=hat.BS.st(ND,fit.JM11,2.5,5)
q7=hat.BS.st(ND,fit.JM11,3,5)
q8=hat.BS.st(ND,fit.JM11,3.5,5)
q9=hat.BS.st(ND,fit.JM11,4,5)
q10=hat.BS.st(ND,fit.JM11,4.5,5)
q11=hat.BS.st(ND,fit.JM11,5,5)
q12=hat.BS.st(ND,fit.JM11,5.5,5)
q13=hat.BS.st(ND,fit.JM11,6,5)
q14=hat.BS.st(ND,fit.JM11,6.5,5)
q15=hat.BS.st(ND,fit.JM11,7,5)
q16=hat.BS.st(ND,fit.JM11,7.5,5)
q17=hat.BS.st(ND,fit.JM11,8,5)
q18=hat.BS.st(ND,fit.JM11,8.5,5)
q19=hat.BS.st(ND,fit.JM11,9,5)
q20=hat.BS.st(ND,fit.JM11,9.5,5)
q21=hat.BS.st(ND,fit.JM11,10,5)



t3.set<-c(z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18,z19,z20,z21)
t5.set<-c(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18,q19,q20,q21)
t3.set
t5.set

plot(y,t5.set,main="blue:t=3 and red:t=5",xlab="s value")
lines(y,t3.set,type="o",lty=2,col="blue")
lines(y,t5.set,pch=22,type="o",lty=3,col="red")
##s =seq(0,10,0.5) with t = 3 or 5



