rm(list=ls())
d=read.csv("CASP.csv")
attach(d)

# install.packages("dplyr")

shapiro.test(RMSD[1:5000])

plot(RMSD,F1)
summary(lm(RMSD~F1+F2+F3+F4+F5+F6+F7+F8+F9))
model=(lm(RMSD~F1+F2+F3+F4+F5+F6+F7+F8+F9))
par(mfrow=c(2,2))
plot(model)
# install.packages("mctest")
library(mctest)
imcdiag(model,method='VIF')
eigprop(model)

#Multicollinearity testing
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




#F1 and F2 is removed

model2=(lm(RMSD~F3.s+F4.s+F5.s+F6.s+F7.s+F8.s+F9.s))
imcdiag(model2,method="VIF")
eigprop(model2)

#F6 is removed

model3=(lm(RMSD~F3.s+F4.s+F5.s+F7.s+F8.s+F9.s))
imcdiag(model3,method="VIF")
eigprop(model3)



X=matrix(c(F1.s,F2.s,F3.s,F4.s,F5.s,F6.s,F7.s,F8.s,F9.s),ncol=9)
X1=(t(X)%*%X)/n
round(X1,3)
eval1=eigen(X1)$values



#eigen value checking for reduced model

X.red=matrix(c(F3.s,F4.s,F5.s,F7.s,F8.s,F9.s),ncol=6)
x1.red=((t(X.red)%*%X.red))/n
eval2=eigen(x1.red)$values



  
#Ridge regression
  
lamda=2.7411
M=X1
for (i in 1:9)
{
  for (j in 1:9)
  {
    if (i==j)
    {
      M[i,j]=M[i,j]+lamda
      
    }
  
      }
  }
  M
  eigen(M)
  

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

par(mfrow=c(1,1))
plot(z,Mat[1,],col='red')
lines(z,Mat[1,],col='red')
points(z,Mat[2,],pch='+',col='blue')
lines(z,Mat[2,],col='blue',lty=8,pch='*')

points(z,Mat[3,],pch='-')
lines(z,Mat[3,],col='green',lty=8,pch='*')

points(z,Mat[4,],pch='-')
lines(z,Mat[4,],col='deep red',lty=8)

points(z,Mat[5,],pch='-')
lines(z,Mat[5,],col='purple',lty=8)

points(z,Mat[6,],pch='-')
lines(z,Mat[6,],col='orange')

points(z,Mat[7,],pch='-')
lines(z,Mat[7,],col='black',lty=8)

points(z,Mat[8,],pch='-')
lines(z,Mat[8,],col='skyblue',lty=8)

lines(z,Mat[9,],pch='-')
points(z,Mat[9,],col='violet',lty=8)














