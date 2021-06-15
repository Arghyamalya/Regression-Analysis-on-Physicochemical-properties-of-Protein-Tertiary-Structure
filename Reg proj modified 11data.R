rm(list=ls())
d=read.csv("CASP2.csv")
attach(d)
# install.packages("dplyr")
shapiro.test(RMSD[1:5000])
length(RMSD)
summary(lm(RMSD~F1+F2+F3+F4+F5+F6+F7+F8+F9))
model=(lm(RMSD~F1+F2+F3+F4+F5+F6+F7+F8+F9))
par(mfrow=c(2,2))
plot(model)
library(mctest)
imcdiag(model,method='VIF')
eigprop(model)

#Multicollinearity Checking

n=length(RMSD)
f=(n-1)/n

F1.s=(F1-mean(F1))/sqrt(var(F1)*f)
F2.s=(F2-mean(F2))/sqrt(var(F2)*f)
F3.s=(F3-mean(F3))/sqrt(var(F3)*f)
F4.s=(F4-mean(F4))/sqrt(var(F4)*f)
F5.s=(F5-mean(F5))/sqrt(var(F5)*f)
F6.s=(F6-mean(F6))/sqrt(var(F6)*f)
F7.s=(F7-mean(F7))/sqrt(var(F7)*f)
F8.s=(F8-mean(F8))/sqrt(var(F8)*f)
F9.s=(F9-mean(F9))/sqrt(var(F9)*f)



model1=(lm(RMSD~F1.s+F2.s+F3.s+F4.s+F5.s+F6.s+F7.s+F8.s+F9.s))
imcdiag(model1,method="VIF")
eigprop(model1)

#(F1,F5) and (F2,F3) are the subsets
# F1 and F2 is removed

model2=(lm(RMSD~F3.s+F4.s+F5.s+F6.s+F7.s+F8.s+F9.s))
imcdiag(model2,method="VIF")
eigprop(model2)


#(F5,F6) subset forms and F6 is removed


model3=(lm(RMSD~F3.s+F4.s+F5.s+F7.s+F8.s+F9.s))
imcdiag(model3,method="VIF")
eigprop(model3)

#(F4,F5) subset and F5 is removed though cond no is<5

model4=(lm(RMSD~F3.s+F4.s+F7.s+F8.s+F9.s))
imcdiag(model4,method="VIF")
eigprop(model4)

#(F4,F9) and F9 is removed

model5=(lm(RMSD~F3.s+F4.s+F7.s+F8.s))
imcdiag(model5,method="VIF")
eigprop(model5)

#no subset forming also all VIFs <15

#eigen values checking~original model

X=matrix(c(F1.s,F2.s,F3.s,F4.s,F5.s,F6.s,F7.s,F8.s,F9.s),ncol=9)
X1=(t(X)%*%X)/n
round(X1,3)
eval1=eigen(X1)$values

#Reduced model

model5=(lm(RMSD~F3.s+F4.s+F7.s+F8.s))
X.red=matrix(c(F3.s,F4.s,F7.s,F8.s),ncol=4)
X1.red=((t(X.red)%*%X.red))/n

rownames(X1.red)=c("F3","F4","F7","F8")
colnames(X1.red)=c("F3","F4","F7","F8")
View(X1.red)
eval2=eigen(x1.red)$values
eval2
summary(model1)







fl=function(x)
{
  m1=solve(X1+(x*diag(9)))%*%t(X)%*%RMSD    
  m1
} 

z=seq(0,3,0.01)  
length(z)
Mat=matrix(0,nrow=9,ncol=length(z))
for(i in 1:length(z))
{
  Mat[,i]=fl(z[i])
}


lam_fin=0.1335939


#trace plot

par(mfrow=c(1,1))
plot(z,Mat[1,],col='red',xlim=c(0,0.5), ylim=c(0,100000),xaxs="i",yaxs="i",
     main="Ridge Trace Plot",col.main="Black"
     ,xlab='lambda',ylab='parameter estimates',text(0.17, 6e+04,expression(lambda==0.1335)))
lines(z,Mat[1,],col='red',lwd=2)
points(z,Mat[2,],col='blue',pch=19)
lines(z,Mat[2,],col='blue',lwd=2)

points(z,Mat[3,],col='green',pch=19)
lines(z,Mat[3,],col='green',lwd=2)

points(z,Mat[4,],col='yellow',pch=19)
lines(z,Mat[4,],col='yellow',lwd=8)

points(z,Mat[5,],col='purple',pch=19)
lines(z,Mat[5,],col='purple',lwd=2)

points(z,Mat[6,],col='orange',pch=19)
lines(z,Mat[6,],col='orange',lwd=6)

points(z,Mat[7,],col='pink',pch=19)
lines(z,Mat[7,],col='pink',lwd=2)

points(z,Mat[8,],col='skyblue',pch=19)
lines(z,Mat[8,],col='skyblue',lwd=2)

lines(z,Mat[9,],col='violet',pch=19)
points(z,Mat[9,],col='violet',lwd=2)
abline(v=lam_fin,col='black',lwd=2)
legend("topright",legend=c("F1","F2","F3","F4","F5","F6","F7","F8","F9"),
       title="corresponding regressors",col=c("red","blue","green","yellow","purple","orange","pink","skyblue","violet"),lty=1,cex=0.7,lwd=3)



#comparison Ridge regression


betar.h=solve((n*X1)+(lam_fin*diag(9)))%*%t(X)%*%RMSD
eval1
evect1=eigen(X1)$vectors
summary(model)
sig.hat=1.66^2
library(Matrix)
 # install.packages("psych")
library(psych)

tr(solve(t(X)%*%X))
m1=matrix(0,nrow=9,ncol=9)
m2=matrix(0,nrow=9,ncol=9)
for(i in 1:9)
{
m1=m1+((evect1[,i]%*%t(evect1[,i]))/eval1[i])
m2=m2+(eval1[i]*(evect1[,i]%*%t(evect1[,i]))/(eval1[i]+lam_fin)^2)
}
cov.betah=sig.hat*(m1)
cov.bridge=sig.hat*(m2)

# sum of variances under usual estimate
sig.hat*tr(m1)

#sum of variances under ridge estimate
sig.hat*tr(m2)

var.usual_e=diag(cov.betah)
var.ridge_e=diag(cov.bridge)
var.dfm=matrix(c(var.usual_e,var.ridge_e),nrow=9)
rownames(var.dfm)=1:9
colnames(var.dfm)=c("Var of usual parameter estimate", "var of Ridge estimate")
View(var.dfm)

cov.betah=sig.hat*(m1)
cov.bridge=sig.hat*(m2)

cov.bridge
cov.betah
tr(m2*sig.hat)+(lam_fin^2)*t(betar.h)%*%solve((n*X1)+(lam_fin*diag(9)))%*%betar.h
mse.ridge=tr(m2*sig.hat)+(lam_fin^2)*t(betar.h)%*%solve((n*X1)+(lam_fin*diag(9)))%*%betar.h
cov.betah
mse.usual=tr(cov.betah)
mse.usual
mse.ridge



da=rnorm(100)
kurtosis((rnorm(500)))
library(e1071)
