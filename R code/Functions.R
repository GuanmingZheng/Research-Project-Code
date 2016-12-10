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




############Draw plot  for JM11

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


t31.set<-c(z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17)
t51.set<-c(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15)
t31.set
t51.set
plot(S,t51.set,main="red:t=5",xlab="s value")
plot(S,t31.set,main="blue:t=3",xlab="s value")

lines(S,t31.set,type="o",lty=2,col="blue")
lines(S,t51.set,pch=22,type="o",lty=3,col="red")
y
S=seq(0,8,0.5)

############Draw plot  for JM15 different method, same model with JM11


S=seq(0,10,0.5)
S
t1=3
t2=5
a1=hat.BS.st(ND,fit.JM15,0,3)
a2=hat.BS.st(ND,fit.JM15,0.5,3)
a3=hat.BS.st(ND,fit.JM15,1,3)
a4=hat.BS.st(ND,fit.JM15,1.5,3)
a5=hat.BS.st(ND,fit.JM15,2,3)
a6=hat.BS.st(ND,fit.JM15,2.5,3)
a7=hat.BS.st(ND,fit.JM15,3,3)
a8=hat.BS.st(ND,fit.JM15,3.5,3)
a9=hat.BS.st(ND,fit.JM15,4,3)
a10=hat.BS.st(ND,fit.JM15,4.5,3)
a11=hat.BS.st(ND,fit.JM15,5,3)
a12=hat.BS.st(ND,fit.JM15,5.5,3)
a13=hat.BS.st(ND,fit.JM15,6,3)
a14=hat.BS.st(ND,fit.JM15,6.5,3)
a15=hat.BS.st(ND,fit.JM15,7,3)
a16=hat.BS.st(ND,fit.JM15,7.5,3)
a17=hat.BS.st(ND,fit.JM15,8,3)
a18=hat.BS.st(ND,fit.JM15,8.5,3)
a19=hat.BS.st(ND,fit.JM15,9,3)
a20=hat.BS.st(ND,fit.JM15,9.5,3)
a21=hat.BS.st(ND,fit.JM15,10,3)

b1=hat.BS.st(ND,fit.JM15,0,5)
b2=hat.BS.st(ND,fit.JM15,0.5,5)
b3=hat.BS.st(ND,fit.JM15,1,5)
b4=hat.BS.st(ND,fit.JM15,1.5,5)
b5=hat.BS.st(ND,fit.JM15,2,5)
b6=hat.BS.st(ND,fit.JM15,2.5,5)
b7=hat.BS.st(ND,fit.JM15,3,5)
b8=hat.BS.st(ND,fit.JM15,3.5,5)
b9=hat.BS.st(ND,fit.JM15,4,5)
b10=hat.BS.st(ND,fit.JM15,4.5,5)
b11=hat.BS.st(ND,fit.JM15,5,5)
b12=hat.BS.st(ND,fit.JM15,5.5,5)
b13=hat.BS.st(ND,fit.JM15,6,5)
b14=hat.BS.st(ND,fit.JM15,6.5,5)
b15=hat.BS.st(ND,fit.JM15,7,5)
b16=hat.BS.st(ND,fit.JM15,7.5,5)
b17=hat.BS.st(ND,fit.JM15,8,5)
b18=hat.BS.st(ND,fit.JM15,8.5,5)
b19=hat.BS.st(ND,fit.JM15,9,5)
b20=hat.BS.st(ND,fit.JM15,9.5,5)
b21=hat.BS.st(ND,fit.JM15,10,5)



##Draw plot for model JM15 t=3 and t=5 
JM15.t3.set<-c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21)
JM15.t5.set<-c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21)

par(mfrow=c(1,1))
plot(S,JM15.t5.set,main="red:t=5",xlab="s value")
lines(S,JM15.t5.set,pch=22,type="o",lty=3,col="red")
plot(S,JM15.t3.set,main="blue:t=3",xlab="s value")
lines(S,JM15.t3.set,type="o",lty=2,col="blue")

## BEFORE LINE BLOW UP 
JM15.t5.set2<-c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15)
S=seq(0,7,0.5)
plot(S,JM15.t5.set2,main="red:t=5",xlab="s value")
lines(S,JM15.t5.set2,pch=22,type="o",lty=3,col="red")

JM15.t3.set2<-c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19)
S=seq(0,9,0.5)
plot(S,JM15.t3.set2,main="blue:t=3",xlab="s value")
lines(S,JM15.t3.set2,type="o",lty=2,col="blue")






##draw plot JM51 same method different model With JM11

c1=hat.BS.st(ND,fit.JM51,0,3)
c2=hat.BS.st(ND,fit.JM51,0.5,3)
c3=hat.BS.st(ND,fit.JM51,1,3)
c4=hat.BS.st(ND,fit.JM51,1.5,3)
c5=hat.BS.st(ND,fit.JM51,2,3)
c6=hat.BS.st(ND,fit.JM51,2.5,3)
c7=hat.BS.st(ND,fit.JM51,3,3)
c8=hat.BS.st(ND,fit.JM51,3.5,3)
c9=hat.BS.st(ND,fit.JM51,4,3)
c10=hat.BS.st(ND,fit.JM51,4.5,3)
c11=hat.BS.st(ND,fit.JM51,5,3)
c12=hat.BS.st(ND,fit.JM51,5.5,3)
c13=hat.BS.st(ND,fit.JM51,6,3)
c14=hat.BS.st(ND,fit.JM51,6.5,3)
c15=hat.BS.st(ND,fit.JM51,7,3)
c16=hat.BS.st(ND,fit.JM51,7.5,3)
c17=hat.BS.st(ND,fit.JM51,8,3)
c18=hat.BS.st(ND,fit.JM51,8.5,3)
c19=hat.BS.st(ND,fit.JM51,9,3)
c20=hat.BS.st(ND,fit.JM51,9.5,3)
c21=hat.BS.st(ND,fit.JM51,10,3)

d1=hat.BS.st(ND,fit.JM51,0,5)
d2=hat.BS.st(ND,fit.JM51,0.5,5)
d3=hat.BS.st(ND,fit.JM51,1,5)
d4=hat.BS.st(ND,fit.JM51,1.5,5)
d5=hat.BS.st(ND,fit.JM51,2,5)
d6=hat.BS.st(ND,fit.JM51,2.5,5)
d7=hat.BS.st(ND,fit.JM51,3,5)
d8=hat.BS.st(ND,fit.JM51,3.5,5)
d9=hat.BS.st(ND,fit.JM51,4,5)
d10=hat.BS.st(ND,fit.JM51,4.5,5)
d11=hat.BS.st(ND,fit.JM51,5,5)
d12=hat.BS.st(ND,fit.JM51,5.5,5)
d13=hat.BS.st(ND,fit.JM51,6,5)
d14=hat.BS.st(ND,fit.JM51,6.5,5)
d15=hat.BS.st(ND,fit.JM51,7,5)
d16=hat.BS.st(ND,fit.JM51,7.5,5)
d17=hat.BS.st(ND,fit.JM51,8,5)
d18=hat.BS.st(ND,fit.JM51,8.5,5)
d19=hat.BS.st(ND,fit.JM51,9,5)
d20=hat.BS.st(ND,fit.JM51,9.5,5)
d21=hat.BS.st(ND,fit.JM51,10,5)


JM51.t3.set<-c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21)
JM51.t5.set<-c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21)

##Draw plot for model JM15 t=3 and t=5 
JM51.t3.set<-c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21)
JM51.t5.set<-c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21)

par(mfrow=c(1,1))
S=seq(0,10,0.5)
plot(S,JM51.t5.set,main="red:t=5",xlab="s value")
lines(S,JM51.t5.set,pch=22,type="o",lty=3,col="red")
plot(S,JM51.t3.set,main="blue:t=3",xlab="s value")
lines(S,JM51.t3.set,type="o",lty=2,col="blue")

## BEFORE LINE BLOW UP 
JM51.t5.set2<-c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15)
S=seq(0,7,0.5)
plot(S,JM51.t5.set2,main="red:t=5",xlab="s value")
lines(S,JM51.t5.set2,pch=22,type="o",lty=3,col="red")

JM51.t3.set2<-c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19)
S=seq(0,9,0.5)
plot(S,JM51.t3.set2,main="blue:t=3",xlab="s value")
lines(S,JM51.t3.set2,type="o",lty=2,col="blue")




