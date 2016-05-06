#!/bin/bash

## Script to reshape mutliple dadi outputs from the modified dadi version ####
## Data are imported (scp *) from the cluster in a desktop computer
## Contains both INDEP and LINKED SIMULATIONS

libR=/RAD_Lamproies/dadi_LAMPROIES/DADI_RESULT/SCRIPT_ANLAYSE_DADI

#Organise the data 
for z in Indep; do for i in AA BET BRE CEN ODO OIR SAU JAL RIS ; do mkdir -p $z/$i/IM $z/$i/AM $z/$i/SC $z/$i/IM2M $z/$i/AM2M $z/$i/SC2M $z/$i/SI ; done ; done

mv *INDEP.fs* Indep/

cd Indep
for j in IM IM2M AM AM2M SC SC2M SI; do mv AA.*_masked_2015_* AA/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv AA_28_28_indep.fs.$j.*_masked_2015_* AA/$j/ ; done 

for j in IM IM2M AM AM2M SC SC2M SI; do mv BET_*INDEP.fs.$j.*_masked_2015_* BET/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv BRE_*INDEP.fs.$j.*_masked_2015_* BRE/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv CEN_*INDEP.fs.$j.*_masked_2015_* CEN/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv ODO_*INDEP.fs.$j.*_masked_2015_* ODO/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv OIR_*INDEP.fs.$j.*_masked_2015_* OIR/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv RIS_*INDEP.fs.$j.*_masked_2015_* RIS/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv SAU_*INDEP.fs.$j.*_masked_2015_* SAU/$j/ ; done 

#Rename folder



#Reshape the results to make it readable in R
for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M AM AM2M SC SC2M SI; do grep -v "Grid" $i.$j.txt | grep "," | sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done
for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M SC SC2M AM AM2M SI; do grep "AIC" $i.$j.txt >> $i.$j.AIC.txt & grep "Optimized log-likelihood" $i.$j.txt >> $i.$j.LogLik.txt & grep "theta" $i.$j.txt >> $i.$j.theta.txt; done ; done
for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M SC SC2M AM AM2M SI; do paste -d " " $i.$j.AIC.txt $i.$j.LogLik.txt $i.$j.theta.txt $i.$j.Opt.Param.txt >> $i.$j.DADI_1.txt  ; done ; done
for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M SC SC2M AM AM2M SI; do sed -e 's/AIC://g' $i.$j.DADI_1.txt | sed -e 's/theta://g' | sed -e 's/Optimized log-likelihood://g' | sed -e 's/ptimized parameters: array//g' | sed -e 's/\]//g' | sed -e 's/\[//g' | sed -e 's/(\|)//g' | sed -e '/nan/d' >> $i.$j.DADI_2.txt ; done ; done
#for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M SC SC2M AM AM2M ; do wc -l $i.$j.DADI_2.txt > $i.$j.Line ; done ; done
#for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M SC SC2M AM AM2M SI; do for `seq in 1 $i.$j.Line` ; do echo $i.$j ; done ; done ; done
rm *Param* *AIC* *DADI_1.txt *theta* *LogLik.txt

mkdir DADI_SYNTHESE_INDEP 
mv *DADI_2.txt DADI_SYNTHESE_INDEP
cd DADI_SYNTHESE_INDEP
#for i in *.DADI_LIE2.txt ; do mv $i ${i%DADI_LIE2.txt}DADI_unfolded.txt ; done

#Compute mean over all files
$libR/dadi_stat_indep_lie.R 

######
paste model DADI.mean.tmp >> mean.tmp.txt
cat entete mean.tmp.txt >> meanFull_Lie.txt
rm model entete *tmp*

#### Now we can edit the 2 files in LibreOffice 
