#At the end, use the script dadi.sh to ad the "-" in the first two colomns
#FS <- read.table("spectrum.OIR.tmp",header=T)
#FS <- read.table("spectrum_SAU.tmp",header=T)
#FS <- read.table("spectrum_CEN.tmp",header=T)
#FS <- read.table("spectrum_ODO.tmp",header=T)
#FS <- read.table("spectrum_AA.tmp2",header=T)
#FS <- read.table("spectrum_BRE.tmp",header=T)
#FS <- read.table("spectrum_RIS.tmp",header=T)

FS <- read.table("spectrum_BET.tmp",header=T)

attach(FS)

NEWFS <- NULL
for (i in 1:nrow(FS))
	{
	FSI <- NULL
	## Only treats monomorphic markers in Petro##
	if ((FS$PEMA[i]+FS$PEMA.1[i]>=4) & (FS$PEMA[i]==0 | FS$PEMA.1[i]==0))
		{
		## If Allele 1 is the majoritary allele in lampetra ##
		if ((FS$LF[i]+FS$LP[i])/(FS$LF[i]+FS$LP[i]+FS$LF.1[i]+FS$LP.1[i])>=0.5)
			{
			## If Allele 1 matches the reference genome leave $lampetra unchanged ##
			if (as.character(FS$Allele1[i]) == substring(as.character(FS$lampe[i]),1))
				{
				FSI <- FS[i,c(1:5,7:9,11:12)]
				## Edit $Punctatus using Allele2 if Allele2 is fixed in the outgroup ##
				if (FS$PEMA.1[i]!=0)
					{
					SubFSI <- substring(as.character(FS$petro[i]),1)
					SubFSI[2] <- as.character(FS$Allele2[i])
					FSI[2] <- paste(SubFSI[2])
					}
				## Or leave $petromyzon unchanged if Allele1 is fixed in the outgroup ##
				else {FSI <- FSI}
				}
			## If Allele 2 matches the reference genome edit $lampetra ##
			if (as.character(FS$Allele2[i]) == substring(as.character(FS$lampe[i]),1))
				{
				FSI <- FS[i,c(1:5,7:9,11:12)]
				## Edit $lampe using Allele2 ##
				SubFSI <- substring(as.character(FS$lampe[i]),1)
				SubFSI[2] <- as.character(FS$Allele1[i])
				FSI[1] <- paste(SubFSI[1],SubFSI[2],SubFSI[3],sep="")
				## Edit $petro using Allele1 if Allele1 is fixed in the outgroup ##
				if (FS$PEMA[i]!=0)
					{
					SubFSI <- substring(as.character(FS$petro[i]),1)
					SubFSI[2] <- as.character(FS$Allele1[i])
					FSI[2] <- paste(SubFSI[2])
					}
				## Or leave $petro unchanged if Allele2 is fixed in the outgroup ##
				else {FSI <- FSI}
				}
			}
		## If Allele 2 is the majoritary allele in lampe ##
		if ((FS$LF[i]+FS$LP[i])/(FS$LF[i]+FS$LP[i]+FS$LF.1[i]+FS$LP.1[i])<0.5)
			{
			## If Allele 1 matches the reference genome edit $lampe ##
			if (as.character(FS$Allele1[i]) == substring(as.character(FS$lampe[i]),1))
				{
				FSI <- FS[i,c(1:5,7:9,11:12)]
				## Edit $lampe using Allele2 ##
				SubFSI <- substring(as.character(FS$lampe[i]),1)
				SubFSI[2] <- as.character(FS$Allele2[i])
				FSI[1] <- paste(SubFSI[1],SubFSI[2],SubFSI[3],sep="")
				## Edit $petro using Allele2 if Allele2 is fixed in the outgroup ##
				if (FS$PEMA.1[i]!=0)
					{
					SubFSI <- substring(as.character(FS$petro[i]),1)
					SubFSI[2] <- as.character(FS$Allele2[i])
					FSI[2] <- paste(SubFSI[2])
					}
				## Or leave $petro unchanged if Allele1 is fixed in the outgroup ##
				else {FSI <- FSI}
				}
			## If Allele 2 matches the reference genome leave $lampe unchanged ##
			if (as.character(FS$Allele2[i]) == substring(as.character(FS$lampe[i]),1))
				{
				FSI <- FS[i,c(1:5,7:9,11:12)]
				## Edit $petro using Allele1 if Allele1 is fixed in the outgroup ##
				if (FS$PEMA[i]!=0)
					{
					SubFSI <- substring(as.character(FS$petro[i]),1)
					SubFSI[2] <- as.character(FS$Allele1[i])
					FSI[2] <- paste(SubFSI[2])
					}
				## Or leave $petro unchanged if Allele2 is fixed in the outgroup ##
				else {FSI <- FSI}
				}
			}
		}
	NEWFS <- rbind(NEWFS,FSI)
	}

z=dim(NEWFS)[1]
x=as.character("-")
t=rep(x,z)
NEWFS[,1]=paste(t,NEWFS[,1],t ,sep="")
NEWFS[,2]=paste(t,NEWFS[,2],t ,sep="")




write.table(NEWFS,"spectrum_BET.data", sep="\t", row.names=F,quote=F)
#write.table(NEWFS,"spectrum_AA.data", sep="\t", row.names=F,quote=F)
#write.table(NEWFS,"spectrum_BRE.data", sep="\t", row.names=F,quote=F)
#write.table(NEWFS,"spectrum_OIR.data", sep="\t", row.names=F,quote=F)
#write.table(NEWFS,"spectrum_SAU.data", sep="\t", row.names=F,quote=F)
#write.table(NEWFS,"spectrum_RIS.data", sep="\t", row.names=F,quote=F)
#write.table(NEWFS,"spectrum_CEN.data", sep="\t", row.names=F,quote=F)
#write.table(NEWFS,"spectrum_ODO.data", sep="\t", row.names=F,quote=F)
