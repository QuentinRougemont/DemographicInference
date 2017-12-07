#!/usr/bin/env Rscript

argv <- commandArgs(T)

vcf  <- argv[1]
mark <- argv[2]

vcf1 <- read.table(vcf)
mark1 <- read.table(mark)

c=merge(mark1, vcf1, by.x="V1",by.y="V3", sort=F)

write.table(c[,c(2:3,1,4:ncol(c))],"vcf_one_snp.tmp",quote=F,row.names=F,col.names=F, sep="\t")

