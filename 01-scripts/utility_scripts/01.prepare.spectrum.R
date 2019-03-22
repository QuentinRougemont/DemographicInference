#!/usr/bin/R

spectrum=read.table("spectrum.tmp", 
                    col.names=c("pop1","pop2","n_pop1",
                    "Allele1", "freq_pop1",
                    "Allele2","freq_pop1",
                    "n_pop2", "freq_pop2","freq_pop2",
                    "scaffold","Position"))

MAJ_pop1 <- round(spectrum[,3]*spectrum[,5],0)
MIN_pop1 <- round(spectrum[,3]*spectrum[,7],0)
MAJ_pop2 <- round(spectrum[,8]*spectrum[,9],0)
MIN_pop2 <- round(spectrum[,8]*spectrum[,10],0)

data=cbind(as.character(spectrum$pop1),as.character(spectrum$pop2),
           as.character(spectrum$Allele1),MAJ_pop1,MAJ_pop2,
           as.character(spectrum$Allele2),
           MIN_pop1,MIN_pop2,
           as.character(spectrum$scaffold),
           spectrum$Position)

colnames(data)=c("pop1","pop2","Allele1","pop1","pop2","Allele2","pop1","pop2","Gene","Position")

#write.table(data, "spectrum.data", col.names=T,row.names=F, quote=F, sep="\t")
#conformation to dadi requirments:
sfs=data
x=as.character("-")
t=rep(x,dim(sfs)[1])
sfs[,1]=paste(t,sfs[,1],t ,sep="")
sfs[,2]=paste(t,sfs[,2],t ,sep="")
write.table(sfs,"spectrum_NP.data", sep="\t", row.names=F,quote=F)

na_rm <- paste("grep -v 'NaN' spectrum_NP.data > spectrum_NP.data2")
system(na_rm, wait = F, intern = F)
