#!/usr/bin/Rscript


rm(list=ls())
all.files=list.files(pattern="DADI_2.txt")

model=c("AA.AM2M","AA.AM","AA.IM2M","AA.IM","AA.SC2M","AA.SC","OIR.AM2M","OIR.AM","OIR.IM2M","OIR.IM","OIR.SC2M","OIR.SC",
       "RIS.AM2M","RIS.AM","RIS.IM2M","RIS.IM","RIS.SC2M","RIS.SC")

entete=c("Riv.Mod","AIC","log_lik","theta","nu1","nu2","m12","m21","me12","me21","Ts","Tsc/am","P") # if oriented add "O"!
#cursor=1
for (f in 1: length(all.files)){
 data=read.table(all.files[f])
 param=subset(data, data[,1]==min(data[,1]))
 write.table((param),"DADI.min.tmp", col.names=F, quote=F, row.names=F, append=T) #col.names=c("AIC","log_lik","theta","nu1","nu2","m12","m21","me12","me21","Ts","Tsc/am","P")
 }
write.table(model, "model", quote=F,col.names=F,row.names=F)
write.table(t(entete),"entete", quote=F,col.names=F,row.names=F)
