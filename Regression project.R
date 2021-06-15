rm(list=ls())
d=read.csv("CASP.csv")
length(d)
attach(d)
plot(F1,RMSD)
plot(lm(RMSD~F1))

