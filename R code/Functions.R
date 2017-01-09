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




####################################START HERE########################################



##JM11-JM16 Functions 
##JM11
S=seq(0,10,0.5)
S
t1=3
t2=5
a1=hat.BS.st(ND,fit.JM11,0,3)
a2=hat.BS.st(ND,fit.JM11,0.5,3)
a3=hat.BS.st(ND,fit.JM11,1,3)
a4=hat.BS.st(ND,fit.JM11,1.5,3)
a5=hat.BS.st(ND,fit.JM11,2,3)
a6=hat.BS.st(ND,fit.JM11,2.5,3)
a7=hat.BS.st(ND,fit.JM11,3,3)
a8=hat.BS.st(ND,fit.JM11,3.5,3)
a9=hat.BS.st(ND,fit.JM11,4,3)
a10=hat.BS.st(ND,fit.JM11,4.5,3)
a11=hat.BS.st(ND,fit.JM11,5,3)
a12=hat.BS.st(ND,fit.JM11,5.5,3)
a13=hat.BS.st(ND,fit.JM11,6,3)
a14=hat.BS.st(ND,fit.JM11,6.5,3)
a15=hat.BS.st(ND,fit.JM11,7,3)
a16=hat.BS.st(ND,fit.JM11,7.5,3)
a17=hat.BS.st(ND,fit.JM11,8,3)
a18=hat.BS.st(ND,fit.JM11,8.5,3)
a19=hat.BS.st(ND,fit.JM11,9,3)
a20=hat.BS.st(ND,fit.JM11,9.5,3)
a21=hat.BS.st(ND,fit.JM11,10,3)

b1=hat.BS.st(ND,fit.JM11,0,5)
b2=hat.BS.st(ND,fit.JM11,0.5,5)
b3=hat.BS.st(ND,fit.JM11,1,5)
b4=hat.BS.st(ND,fit.JM11,1.5,5)
b5=hat.BS.st(ND,fit.JM11,2,5)
b6=hat.BS.st(ND,fit.JM11,2.5,5)
b7=hat.BS.st(ND,fit.JM11,3,5)
b8=hat.BS.st(ND,fit.JM11,3.5,5)
b9=hat.BS.st(ND,fit.JM11,4,5)
b10=hat.BS.st(ND,fit.JM11,4.5,5)
b11=hat.BS.st(ND,fit.JM11,5,5)
b12=hat.BS.st(ND,fit.JM11,5.5,5)
b13=hat.BS.st(ND,fit.JM11,6,5)
b14=hat.BS.st(ND,fit.JM11,6.5,5)
b15=hat.BS.st(ND,fit.JM11,7,5)
b16=hat.BS.st(ND,fit.JM11,7.5,5)
b17=hat.BS.st(ND,fit.JM11,8,5)
b18=hat.BS.st(ND,fit.JM11,8.5,5)
b19=hat.BS.st(ND,fit.JM11,9,5)
b20=hat.BS.st(ND,fit.JM11,9.5,5)
b21=hat.BS.st(ND,fit.JM11,10,5)



##Draw plot for model JM15 t=3 and t=5 
JM11.t3.set<-c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21)
JM11.t5.set<-c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21)



##JM12
S=seq(0,10,0.5)
S
t1=3
t2=5
c1=hat.BS.st(ND,fit.JM12,0,3)
c2=hat.BS.st(ND,fit.JM12,0.5,3)
c3=hat.BS.st(ND,fit.JM12,1,3)
c4=hat.BS.st(ND,fit.JM12,1.5,3)
c5=hat.BS.st(ND,fit.JM12,2,3)
c6=hat.BS.st(ND,fit.JM12,2.5,3)
c7=hat.BS.st(ND,fit.JM12,3,3)
c8=hat.BS.st(ND,fit.JM12,3.5,3)
c9=hat.BS.st(ND,fit.JM12,4,3)
c10=hat.BS.st(ND,fit.JM12,4.5,3)
c11=hat.BS.st(ND,fit.JM12,5,3)
c12=hat.BS.st(ND,fit.JM12,5.5,3)
c13=hat.BS.st(ND,fit.JM12,6,3)
c14=hat.BS.st(ND,fit.JM12,6.5,3)
c15=hat.BS.st(ND,fit.JM12,7,3)
c16=hat.BS.st(ND,fit.JM12,7.5,3)
c17=hat.BS.st(ND,fit.JM12,8,3)
c18=hat.BS.st(ND,fit.JM12,8.5,3)
c19=hat.BS.st(ND,fit.JM12,9,3)
c20=hat.BS.st(ND,fit.JM12,9.5,3)
c21=hat.BS.st(ND,fit.JM12,10,3)

d1=hat.BS.st(ND,fit.JM12,0,5)
d2=hat.BS.st(ND,fit.JM12,0.5,5)
d3=hat.BS.st(ND,fit.JM12,1,5)
d4=hat.BS.st(ND,fit.JM12,1.5,5)
d5=hat.BS.st(ND,fit.JM12,2,5)
d6=hat.BS.st(ND,fit.JM12,2.5,5)
d7=hat.BS.st(ND,fit.JM12,3,5)
d8=hat.BS.st(ND,fit.JM12,3.5,5)
d9=hat.BS.st(ND,fit.JM12,4,5)
d10=hat.BS.st(ND,fit.JM12,4.5,5)
d11=hat.BS.st(ND,fit.JM12,5,5)
d12=hat.BS.st(ND,fit.JM12,5.5,5)
d13=hat.BS.st(ND,fit.JM12,6,5)
d14=hat.BS.st(ND,fit.JM12,6.5,5)
d15=hat.BS.st(ND,fit.JM12,7,5)
d16=hat.BS.st(ND,fit.JM12,7.5,5)
d17=hat.BS.st(ND,fit.JM12,8,5)
d18=hat.BS.st(ND,fit.JM12,8.5,5)
d19=hat.BS.st(ND,fit.JM12,9,5)
d20=hat.BS.st(ND,fit.JM12,9.5,5)
d21=hat.BS.st(ND,fit.JM12,10,5)


##Draw plot for model JM15 t=3 and t=5 
JM12.t3.set<-c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21)
JM12.t5.set<-c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21)


#JM13

S=seq(0,10,0.5)
S
t1=3
t2=5
e1=hat.BS.st(ND,fit.JM13,0,3)
e2=hat.BS.st(ND,fit.JM13,0.5,3)
e3=hat.BS.st(ND,fit.JM13,1,3)
e4=hat.BS.st(ND,fit.JM13,1.5,3)
e5=hat.BS.st(ND,fit.JM13,2,3)
e6=hat.BS.st(ND,fit.JM13,2.5,3)
e7=hat.BS.st(ND,fit.JM13,3,3)
e8=hat.BS.st(ND,fit.JM13,3.5,3)
e9=hat.BS.st(ND,fit.JM13,4,3)
e10=hat.BS.st(ND,fit.JM13,4.5,3)
e11=hat.BS.st(ND,fit.JM13,5,3)
e12=hat.BS.st(ND,fit.JM13,5.5,3)
e13=hat.BS.st(ND,fit.JM13,6,3)
e14=hat.BS.st(ND,fit.JM13,6.5,3)
e15=hat.BS.st(ND,fit.JM13,7,3)
e16=hat.BS.st(ND,fit.JM13,7.5,3)
e17=hat.BS.st(ND,fit.JM13,8,3)
e18=hat.BS.st(ND,fit.JM13,8.5,3)
e19=hat.BS.st(ND,fit.JM13,9,3)
e20=hat.BS.st(ND,fit.JM13,9.5,3)
e21=hat.BS.st(ND,fit.JM13,10,3)

f1=hat.BS.st(ND,fit.JM13,0,5)
f2=hat.BS.st(ND,fit.JM13,0.5,5)
f3=hat.BS.st(ND,fit.JM13,1,5)
f4=hat.BS.st(ND,fit.JM13,1.5,5)
f5=hat.BS.st(ND,fit.JM13,2,5)
f6=hat.BS.st(ND,fit.JM13,2.5,5)
f7=hat.BS.st(ND,fit.JM13,3,5)
f8=hat.BS.st(ND,fit.JM13,3.5,5)
f9=hat.BS.st(ND,fit.JM13,4,5)
f10=hat.BS.st(ND,fit.JM13,4.5,5)
f11=hat.BS.st(ND,fit.JM13,5,5)
f12=hat.BS.st(ND,fit.JM13,5.5,5)
f13=hat.BS.st(ND,fit.JM13,6,5)
f14=hat.BS.st(ND,fit.JM13,6.5,5)
f15=hat.BS.st(ND,fit.JM13,7,5)
f16=hat.BS.st(ND,fit.JM13,7.5,5)
f17=hat.BS.st(ND,fit.JM13,8,5)
f18=hat.BS.st(ND,fit.JM13,8.5,5)
f19=hat.BS.st(ND,fit.JM13,9,5)
f20=hat.BS.st(ND,fit.JM13,9.5,5)
f21=hat.BS.st(ND,fit.JM13,10,5)


JM13.t3.set<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20,e21)
JM13.t5.set<-c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21)


#JM14
S=seq(0,10,0.5)
S
t1=3
t2=5
g1=hat.BS.st(ND,fit.JM14,0,3)
g2=hat.BS.st(ND,fit.JM14,0.5,3)
g3=hat.BS.st(ND,fit.JM14,1,3)
g4=hat.BS.st(ND,fit.JM14,1.5,3)
g5=hat.BS.st(ND,fit.JM14,2,3)
g6=hat.BS.st(ND,fit.JM14,2.5,3)
g7=hat.BS.st(ND,fit.JM14,3,3)
g8=hat.BS.st(ND,fit.JM14,3.5,3)
g9=hat.BS.st(ND,fit.JM14,4,3)
g10=hat.BS.st(ND,fit.JM14,4.5,3)
g11=hat.BS.st(ND,fit.JM14,5,3)
g12=hat.BS.st(ND,fit.JM14,5.5,3)
g13=hat.BS.st(ND,fit.JM14,6,3)
g14=hat.BS.st(ND,fit.JM14,6.5,3)
g15=hat.BS.st(ND,fit.JM14,7,3)
g16=hat.BS.st(ND,fit.JM14,7.5,3)
g17=hat.BS.st(ND,fit.JM14,8,3)
g18=hat.BS.st(ND,fit.JM14,8.5,3)
g19=hat.BS.st(ND,fit.JM14,9,3)
g20=hat.BS.st(ND,fit.JM14,9.5,3)
g21=hat.BS.st(ND,fit.JM14,10,3)

h1=hat.BS.st(ND,fit.JM14,0,5)
h2=hat.BS.st(ND,fit.JM14,0.5,5)
h3=hat.BS.st(ND,fit.JM14,1,5)
h4=hat.BS.st(ND,fit.JM14,1.5,5)
h5=hat.BS.st(ND,fit.JM14,2,5)
h6=hat.BS.st(ND,fit.JM14,2.5,5)
h7=hat.BS.st(ND,fit.JM14,3,5)
h8=hat.BS.st(ND,fit.JM14,3.5,5)
h9=hat.BS.st(ND,fit.JM14,4,5)
h10=hat.BS.st(ND,fit.JM14,4.5,5)
h11=hat.BS.st(ND,fit.JM14,5,5)
h12=hat.BS.st(ND,fit.JM14,5.5,5)
h13=hat.BS.st(ND,fit.JM14,6,5)
h14=hat.BS.st(ND,fit.JM14,6.5,5)
h15=hat.BS.st(ND,fit.JM14,7,5)
h16=hat.BS.st(ND,fit.JM14,7.5,5)
h17=hat.BS.st(ND,fit.JM14,8,5)
h18=hat.BS.st(ND,fit.JM14,8.5,5)
h19=hat.BS.st(ND,fit.JM14,9,5)
h20=hat.BS.st(ND,fit.JM14,9.5,5)
h21=hat.BS.st(ND,fit.JM14,10,5)


JM14.t3.set<-c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21)
JM14.t5.set<-c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19,h20,h21)


##JM15

S=seq(0,10,0.5)
S
t1=3
t2=5
i1=hat.BS.st(ND,fit.JM15,0,3)
i2=hat.BS.st(ND,fit.JM15,0.5,3)
i3=hat.BS.st(ND,fit.JM15,1,3)
i4=hat.BS.st(ND,fit.JM15,1.5,3)
i5=hat.BS.st(ND,fit.JM15,2,3)
i6=hat.BS.st(ND,fit.JM15,2.5,3)
i7=hat.BS.st(ND,fit.JM15,3,3)
i8=hat.BS.st(ND,fit.JM15,3.5,3)
i9=hat.BS.st(ND,fit.JM15,4,3)
i10=hat.BS.st(ND,fit.JM15,4.5,3)
i11=hat.BS.st(ND,fit.JM15,5,3)
i12=hat.BS.st(ND,fit.JM15,5.5,3)
i13=hat.BS.st(ND,fit.JM15,6,3)
i14=hat.BS.st(ND,fit.JM15,6.5,3)
i15=hat.BS.st(ND,fit.JM15,7,3)
i16=hat.BS.st(ND,fit.JM15,7.5,3)
i17=hat.BS.st(ND,fit.JM15,8,3)
i18=hat.BS.st(ND,fit.JM15,8.5,3)
i19=hat.BS.st(ND,fit.JM15,9,3)
i20=hat.BS.st(ND,fit.JM15,9.5,3)
i21=hat.BS.st(ND,fit.JM15,10,3)

j1=hat.BS.st(ND,fit.JM15,0,5)
j2=hat.BS.st(ND,fit.JM15,0.5,5)
j3=hat.BS.st(ND,fit.JM15,1,5)
j4=hat.BS.st(ND,fit.JM15,1.5,5)
j5=hat.BS.st(ND,fit.JM15,2,5)
j6=hat.BS.st(ND,fit.JM15,2.5,5)
j7=hat.BS.st(ND,fit.JM15,3,5)
j8=hat.BS.st(ND,fit.JM15,3.5,5)
j9=hat.BS.st(ND,fit.JM15,4,5)
j10=hat.BS.st(ND,fit.JM15,4.5,5)
j11=hat.BS.st(ND,fit.JM15,5,5)
j12=hat.BS.st(ND,fit.JM15,5.5,5)
j13=hat.BS.st(ND,fit.JM15,6,5)
j14=hat.BS.st(ND,fit.JM15,6.5,5)
j15=hat.BS.st(ND,fit.JM15,7,5)
j16=hat.BS.st(ND,fit.JM15,7.5,5)
j17=hat.BS.st(ND,fit.JM15,8,5)
j18=hat.BS.st(ND,fit.JM15,8.5,5)
j19=hat.BS.st(ND,fit.JM15,9,5)
j20=hat.BS.st(ND,fit.JM15,9.5,5)
j21=hat.BS.st(ND,fit.JM15,10,5)


JM15.t3.set<-c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21)
JM15.t5.set<-c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,j17,j18,j19,j20,j21)



##JM21-JM26


##JM21
S=seq(0,10,0.5)
S
t1=3
t2=5
k1=hat.BS.st(ND,fit.JM21,0,3)
k2=hat.BS.st(ND,fit.JM21,0.5,3)
k3=hat.BS.st(ND,fit.JM21,1,3)
k4=hat.BS.st(ND,fit.JM21,1.5,3)
k5=hat.BS.st(ND,fit.JM21,2,3)
k6=hat.BS.st(ND,fit.JM21,2.5,3)
k7=hat.BS.st(ND,fit.JM21,3,3)
k8=hat.BS.st(ND,fit.JM21,3.5,3)
k9=hat.BS.st(ND,fit.JM21,4,3)
k10=hat.BS.st(ND,fit.JM21,4.5,3)
k11=hat.BS.st(ND,fit.JM21,5,3)
k12=hat.BS.st(ND,fit.JM21,5.5,3)
k13=hat.BS.st(ND,fit.JM21,6,3)
k14=hat.BS.st(ND,fit.JM21,6.5,3)
k15=hat.BS.st(ND,fit.JM21,7,3)
k16=hat.BS.st(ND,fit.JM21,7.5,3)
k17=hat.BS.st(ND,fit.JM21,8,3)
k18=hat.BS.st(ND,fit.JM21,8.5,3)
k19=hat.BS.st(ND,fit.JM21,9,3)
k20=hat.BS.st(ND,fit.JM21,9.5,3)
k21=hat.BS.st(ND,fit.JM21,10,3)

l1=hat.BS.st(ND,fit.JM21,0,5)
l2=hat.BS.st(ND,fit.JM21,0.5,5)
l3=hat.BS.st(ND,fit.JM21,1,5)
l4=hat.BS.st(ND,fit.JM21,1.5,5)
l5=hat.BS.st(ND,fit.JM21,2,5)
l6=hat.BS.st(ND,fit.JM21,2.5,5)
l7=hat.BS.st(ND,fit.JM21,3,5)
l8=hat.BS.st(ND,fit.JM21,3.5,5)
l9=hat.BS.st(ND,fit.JM21,4,5)
l10=hat.BS.st(ND,fit.JM21,4.5,5)
l11=hat.BS.st(ND,fit.JM21,5,5)
l12=hat.BS.st(ND,fit.JM21,5.5,5)
l13=hat.BS.st(ND,fit.JM21,6,5)
l14=hat.BS.st(ND,fit.JM21,6.5,5)
l15=hat.BS.st(ND,fit.JM21,7,5)
l16=hat.BS.st(ND,fit.JM21,7.5,5)
l17=hat.BS.st(ND,fit.JM21,8,5)
l18=hat.BS.st(ND,fit.JM21,8.5,5)
l19=hat.BS.st(ND,fit.JM21,9,5)
l20=hat.BS.st(ND,fit.JM21,9.5,5)
l21=hat.BS.st(ND,fit.JM21,10,5)


JM21.t3.set<-c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21)
JM21.t5.set<-c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21)

##JM22

S=seq(0,10,0.5)
S
t1=3
t2=5
m1=hat.BS.st(ND,fit.JM22,0,3)
m2=hat.BS.st(ND,fit.JM22,0.5,3)
m3=hat.BS.st(ND,fit.JM22,1,3)
m4=hat.BS.st(ND,fit.JM22,1.5,3)
m5=hat.BS.st(ND,fit.JM22,2,3)
m6=hat.BS.st(ND,fit.JM22,2.5,3)
m7=hat.BS.st(ND,fit.JM22,3,3)
m8=hat.BS.st(ND,fit.JM22,3.5,3)
m9=hat.BS.st(ND,fit.JM22,4,3)
m10=hat.BS.st(ND,fit.JM22,4.5,3)
m11=hat.BS.st(ND,fit.JM22,5,3)
m12=hat.BS.st(ND,fit.JM22,5.5,3)
m13=hat.BS.st(ND,fit.JM22,6,3)
m14=hat.BS.st(ND,fit.JM22,6.5,3)
m15=hat.BS.st(ND,fit.JM22,7,3)
m16=hat.BS.st(ND,fit.JM22,7.5,3)
m17=hat.BS.st(ND,fit.JM22,8,3)
m18=hat.BS.st(ND,fit.JM22,8.5,3)
m19=hat.BS.st(ND,fit.JM22,9,3)
m20=hat.BS.st(ND,fit.JM22,9.5,3)
m21=hat.BS.st(ND,fit.JM22,10,3)

n1=hat.BS.st(ND,fit.JM22,0,5)
n2=hat.BS.st(ND,fit.JM22,0.5,5)
n3=hat.BS.st(ND,fit.JM22,1,5)
n4=hat.BS.st(ND,fit.JM22,1.5,5)
n5=hat.BS.st(ND,fit.JM22,2,5)
n6=hat.BS.st(ND,fit.JM22,2.5,5)
n7=hat.BS.st(ND,fit.JM22,3,5)
n8=hat.BS.st(ND,fit.JM22,3.5,5)
n9=hat.BS.st(ND,fit.JM22,4,5)
n10=hat.BS.st(ND,fit.JM22,4.5,5)
n11=hat.BS.st(ND,fit.JM22,5,5)
n12=hat.BS.st(ND,fit.JM22,5.5,5)
n13=hat.BS.st(ND,fit.JM22,6,5)
n14=hat.BS.st(ND,fit.JM22,6.5,5)
n15=hat.BS.st(ND,fit.JM22,7,5)
n16=hat.BS.st(ND,fit.JM22,7.5,5)
n17=hat.BS.st(ND,fit.JM22,8,5)
n18=hat.BS.st(ND,fit.JM22,8.5,5)
n19=hat.BS.st(ND,fit.JM22,9,5)
n20=hat.BS.st(ND,fit.JM22,9.5,5)
n21=hat.BS.st(ND,fit.JM22,10,5)


JM22.t3.set<-c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21)
JM22.t5.set<-c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,n21)


##JM23
S=seq(0,10,0.5)
S
t1=3
t2=5
o1=hat.BS.st(ND,fit.JM23,0,3)
o2=hat.BS.st(ND,fit.JM23,0.5,3)
o3=hat.BS.st(ND,fit.JM23,1,3)
o4=hat.BS.st(ND,fit.JM23,1.5,3)
o5=hat.BS.st(ND,fit.JM23,2,3)
o6=hat.BS.st(ND,fit.JM23,2.5,3)
o7=hat.BS.st(ND,fit.JM23,3,3)
o8=hat.BS.st(ND,fit.JM23,3.5,3)
o9=hat.BS.st(ND,fit.JM23,4,3)
o10=hat.BS.st(ND,fit.JM23,4.5,3)
o11=hat.BS.st(ND,fit.JM23,5,3)
o12=hat.BS.st(ND,fit.JM23,5.5,3)
o13=hat.BS.st(ND,fit.JM23,6,3)
o14=hat.BS.st(ND,fit.JM23,6.5,3)
o15=hat.BS.st(ND,fit.JM23,7,3)
o16=hat.BS.st(ND,fit.JM23,7.5,3)
o17=hat.BS.st(ND,fit.JM23,8,3)
o18=hat.BS.st(ND,fit.JM23,8.5,3)
o19=hat.BS.st(ND,fit.JM23,9,3)
o20=hat.BS.st(ND,fit.JM23,9.5,3)
o21=hat.BS.st(ND,fit.JM23,10,3)

p1=hat.BS.st(ND,fit.JM23,0,5)
p2=hat.BS.st(ND,fit.JM23,0.5,5)
p3=hat.BS.st(ND,fit.JM23,1,5)
p4=hat.BS.st(ND,fit.JM23,1.5,5)
p5=hat.BS.st(ND,fit.JM23,2,5)
p6=hat.BS.st(ND,fit.JM23,2.5,5)
p7=hat.BS.st(ND,fit.JM23,3,5)
p8=hat.BS.st(ND,fit.JM23,3.5,5)
p9=hat.BS.st(ND,fit.JM23,4,5)
p10=hat.BS.st(ND,fit.JM23,4.5,5)
p11=hat.BS.st(ND,fit.JM23,5,5)
p12=hat.BS.st(ND,fit.JM23,5.5,5)
p13=hat.BS.st(ND,fit.JM23,6,5)
p14=hat.BS.st(ND,fit.JM23,6.5,5)
p15=hat.BS.st(ND,fit.JM23,7,5)
p16=hat.BS.st(ND,fit.JM23,7.5,5)
p17=hat.BS.st(ND,fit.JM23,8,5)
p18=hat.BS.st(ND,fit.JM23,8.5,5)
p19=hat.BS.st(ND,fit.JM23,9,5)
p20=hat.BS.st(ND,fit.JM23,9.5,5)
p21=hat.BS.st(ND,fit.JM23,10,5)


JM23.t3.set<-c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13,o14,o15,o16,o17,o18,o19,o20,o21)
JM23.t5.set<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21)


##JM24
S=seq(0,10,0.5)
S
t1=3
t2=5
q1=hat.BS.st(ND,fit.JM24,0,3)
q2=hat.BS.st(ND,fit.JM24,0.5,3)
q3=hat.BS.st(ND,fit.JM24,1,3)
q4=hat.BS.st(ND,fit.JM24,1.5,3)
q5=hat.BS.st(ND,fit.JM24,2,3)
q6=hat.BS.st(ND,fit.JM24,2.5,3)
q7=hat.BS.st(ND,fit.JM24,3,3)
q8=hat.BS.st(ND,fit.JM24,3.5,3)
q9=hat.BS.st(ND,fit.JM24,4,3)
q10=hat.BS.st(ND,fit.JM24,4.5,3)
q11=hat.BS.st(ND,fit.JM24,5,3)
q12=hat.BS.st(ND,fit.JM24,5.5,3)
q13=hat.BS.st(ND,fit.JM24,6,3)
q14=hat.BS.st(ND,fit.JM24,6.5,3)
q15=hat.BS.st(ND,fit.JM24,7,3)
q16=hat.BS.st(ND,fit.JM24,7.5,3)
q17=hat.BS.st(ND,fit.JM24,8,3)
q18=hat.BS.st(ND,fit.JM24,8.5,3)
q19=hat.BS.st(ND,fit.JM24,9,3)
q20=hat.BS.st(ND,fit.JM24,9.5,3)
q21=hat.BS.st(ND,fit.JM24,10,3)

r1=hat.BS.st(ND,fit.JM24,0,5)
r2=hat.BS.st(ND,fit.JM24,0.5,5)
r3=hat.BS.st(ND,fit.JM24,1,5)
r4=hat.BS.st(ND,fit.JM24,1.5,5)
r5=hat.BS.st(ND,fit.JM24,2,5)
r6=hat.BS.st(ND,fit.JM24,2.5,5)
r7=hat.BS.st(ND,fit.JM24,3,5)
r8=hat.BS.st(ND,fit.JM24,3.5,5)
r9=hat.BS.st(ND,fit.JM24,4,5)
r10=hat.BS.st(ND,fit.JM24,4.5,5)
r11=hat.BS.st(ND,fit.JM24,5,5)
r12=hat.BS.st(ND,fit.JM24,5.5,5)
r13=hat.BS.st(ND,fit.JM24,6,5)
r14=hat.BS.st(ND,fit.JM24,6.5,5)
r15=hat.BS.st(ND,fit.JM24,7,5)
r16=hat.BS.st(ND,fit.JM24,7.5,5)
r17=hat.BS.st(ND,fit.JM24,8,5)
r18=hat.BS.st(ND,fit.JM24,8.5,5)
r19=hat.BS.st(ND,fit.JM24,9,5)
r20=hat.BS.st(ND,fit.JM24,9.5,5)
r21=hat.BS.st(ND,fit.JM24,10,5)


JM24.t3.set<-c(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18,q19,q20,q21)
JM24.t5.set<-c(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21)



##JM25
S=seq(0,10,0.5)
S
t1=3
t2=5
s1=hat.BS.st(ND,fit.JM25,0,3)
s2=hat.BS.st(ND,fit.JM25,0.5,3)
s3=hat.BS.st(ND,fit.JM25,1,3)
s4=hat.BS.st(ND,fit.JM25,1.5,3)
s5=hat.BS.st(ND,fit.JM25,2,3)
s6=hat.BS.st(ND,fit.JM25,2.5,3)
s7=hat.BS.st(ND,fit.JM25,3,3)
s8=hat.BS.st(ND,fit.JM25,3.5,3)
s9=hat.BS.st(ND,fit.JM25,4,3)
s10=hat.BS.st(ND,fit.JM25,4.5,3)
s11=hat.BS.st(ND,fit.JM25,5,3)
s12=hat.BS.st(ND,fit.JM25,5.5,3)
s13=hat.BS.st(ND,fit.JM25,6,3)
s14=hat.BS.st(ND,fit.JM25,6.5,3)
s15=hat.BS.st(ND,fit.JM25,7,3)
s16=hat.BS.st(ND,fit.JM25,7.5,3)
s17=hat.BS.st(ND,fit.JM25,8,3)
s18=hat.BS.st(ND,fit.JM25,8.5,3)
s19=hat.BS.st(ND,fit.JM25,9,3)
s20=hat.BS.st(ND,fit.JM25,9.5,3)
s21=hat.BS.st(ND,fit.JM25,10,3)

t1=hat.BS.st(ND,fit.JM25,0,5)
t2=hat.BS.st(ND,fit.JM25,0.5,5)
t3=hat.BS.st(ND,fit.JM25,1,5)
t4=hat.BS.st(ND,fit.JM25,1.5,5)
t5=hat.BS.st(ND,fit.JM25,2,5)
t6=hat.BS.st(ND,fit.JM25,2.5,5)
t7=hat.BS.st(ND,fit.JM25,3,5)
t8=hat.BS.st(ND,fit.JM25,3.5,5)
t9=hat.BS.st(ND,fit.JM25,4,5)
t10=hat.BS.st(ND,fit.JM25,4.5,5)
t11=hat.BS.st(ND,fit.JM25,5,5)
t12=hat.BS.st(ND,fit.JM25,5.5,5)
t13=hat.BS.st(ND,fit.JM25,6,5)
t14=hat.BS.st(ND,fit.JM25,6.5,5)
t15=hat.BS.st(ND,fit.JM25,7,5)
t16=hat.BS.st(ND,fit.JM25,7.5,5)
t17=hat.BS.st(ND,fit.JM25,8,5)
t18=hat.BS.st(ND,fit.JM25,8.5,5)
t19=hat.BS.st(ND,fit.JM25,9,5)
t20=hat.BS.st(ND,fit.JM25,9.5,5)
t21=hat.BS.st(ND,fit.JM25,10,5)


JM25.t3.set<-c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21)
JM25.t5.set<-c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21)


##JM31-35
##JM31
S=seq(0,10,0.5)
S
t1=3
t2=5
u1=hat.BS.st(ND,fit.JM31,0,3)
u2=hat.BS.st(ND,fit.JM31,0.5,3)
u3=hat.BS.st(ND,fit.JM31,1,3)
u4=hat.BS.st(ND,fit.JM31,1.5,3)
u5=hat.BS.st(ND,fit.JM31,2,3)
u6=hat.BS.st(ND,fit.JM31,2.5,3)
u7=hat.BS.st(ND,fit.JM31,3,3)
u8=hat.BS.st(ND,fit.JM31,3.5,3)
u9=hat.BS.st(ND,fit.JM31,4,3)
u10=hat.BS.st(ND,fit.JM31,4.5,3)
u11=hat.BS.st(ND,fit.JM31,5,3)
u12=hat.BS.st(ND,fit.JM31,5.5,3)
u13=hat.BS.st(ND,fit.JM31,6,3)
u14=hat.BS.st(ND,fit.JM31,6.5,3)
u15=hat.BS.st(ND,fit.JM31,7,3)
u16=hat.BS.st(ND,fit.JM31,7.5,3)
u17=hat.BS.st(ND,fit.JM31,8,3)
u18=hat.BS.st(ND,fit.JM31,8.5,3)
u19=hat.BS.st(ND,fit.JM31,9,3)
u20=hat.BS.st(ND,fit.JM31,9.5,3)
u21=hat.BS.st(ND,fit.JM31,10,3)

v1=hat.BS.st(ND,fit.JM31,0,5)
v2=hat.BS.st(ND,fit.JM31,0.5,5)
v3=hat.BS.st(ND,fit.JM31,1,5)
v4=hat.BS.st(ND,fit.JM31,1.5,5)
v5=hat.BS.st(ND,fit.JM31,2,5)
v6=hat.BS.st(ND,fit.JM31,2.5,5)
v7=hat.BS.st(ND,fit.JM31,3,5)
v8=hat.BS.st(ND,fit.JM31,3.5,5)
v9=hat.BS.st(ND,fit.JM31,4,5)
v10=hat.BS.st(ND,fit.JM31,4.5,5)
v11=hat.BS.st(ND,fit.JM31,5,5)
v12=hat.BS.st(ND,fit.JM31,5.5,5)
v13=hat.BS.st(ND,fit.JM31,6,5)
v14=hat.BS.st(ND,fit.JM31,6.5,5)
v15=hat.BS.st(ND,fit.JM31,7,5)
v16=hat.BS.st(ND,fit.JM31,7.5,5)
v17=hat.BS.st(ND,fit.JM31,8,5)
v18=hat.BS.st(ND,fit.JM31,8.5,5)
v19=hat.BS.st(ND,fit.JM31,9,5)
v20=hat.BS.st(ND,fit.JM31,9.5,5)
v21=hat.BS.st(ND,fit.JM31,10,5)


JM31.t3.set<-c(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21)
JM31.t5.set<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20,v21)


#JM32

S=seq(0,10,0.5)
S
t1=3
t2=5
u1=hat.BS.st(ND,fit.JM31,0,3)
u2=hat.BS.st(ND,fit.JM31,0.5,3)
u3=hat.BS.st(ND,fit.JM31,1,3)
u4=hat.BS.st(ND,fit.JM31,1.5,3)
u5=hat.BS.st(ND,fit.JM31,2,3)
u6=hat.BS.st(ND,fit.JM31,2.5,3)
u7=hat.BS.st(ND,fit.JM31,3,3)
u8=hat.BS.st(ND,fit.JM31,3.5,3)
u9=hat.BS.st(ND,fit.JM31,4,3)
u10=hat.BS.st(ND,fit.JM31,4.5,3)
u11=hat.BS.st(ND,fit.JM31,5,3)
u12=hat.BS.st(ND,fit.JM31,5.5,3)
u13=hat.BS.st(ND,fit.JM31,6,3)
u14=hat.BS.st(ND,fit.JM31,6.5,3)
u15=hat.BS.st(ND,fit.JM31,7,3)
u16=hat.BS.st(ND,fit.JM31,7.5,3)
u17=hat.BS.st(ND,fit.JM31,8,3)
u18=hat.BS.st(ND,fit.JM31,8.5,3)
u19=hat.BS.st(ND,fit.JM31,9,3)
u20=hat.BS.st(ND,fit.JM31,9.5,3)
u21=hat.BS.st(ND,fit.JM31,10,3)

v1=hat.BS.st(ND,fit.JM31,0,5)
v2=hat.BS.st(ND,fit.JM31,0.5,5)
v3=hat.BS.st(ND,fit.JM31,1,5)
v4=hat.BS.st(ND,fit.JM31,1.5,5)
v5=hat.BS.st(ND,fit.JM31,2,5)
v6=hat.BS.st(ND,fit.JM31,2.5,5)
v7=hat.BS.st(ND,fit.JM31,3,5)
v8=hat.BS.st(ND,fit.JM31,3.5,5)
v9=hat.BS.st(ND,fit.JM31,4,5)
v10=hat.BS.st(ND,fit.JM31,4.5,5)
v11=hat.BS.st(ND,fit.JM31,5,5)
v12=hat.BS.st(ND,fit.JM31,5.5,5)
v13=hat.BS.st(ND,fit.JM31,6,5)
v14=hat.BS.st(ND,fit.JM31,6.5,5)
v15=hat.BS.st(ND,fit.JM31,7,5)
v16=hat.BS.st(ND,fit.JM31,7.5,5)
v17=hat.BS.st(ND,fit.JM31,8,5)
v18=hat.BS.st(ND,fit.JM31,8.5,5)
v19=hat.BS.st(ND,fit.JM31,9,5)
v20=hat.BS.st(ND,fit.JM31,9.5,5)
v21=hat.BS.st(ND,fit.JM31,10,5)


JM31.t3.set<-c(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21)
JM31.t5.set<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20,v21)


