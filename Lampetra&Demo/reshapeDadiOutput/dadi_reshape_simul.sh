#!/bin/bash

## Script to reshape mutliple dadi outputs from the modified dadi version ####
## Data are imported (scp *) from the cluster in a desktop computer
## Contains both INDEP and LINKED SIMULATIONS

libR=/RAD_Lamproies/dadi_LAMPROIES/DADI_RESULT/SCRIPT_ANLAYSE_DADI

#Organise the data 
for z in Indep; do for i in AA BET BRE CEN ODO OIR SAU JAL RIS ; do mkdir -p $z/$i/IM $z/$i/AM $z/$i/SC $z/$i/IM2M $z/$i/AM2M $z/$i/SC2M $z/$i/SI ; done ; done
for z in PEMA ; do for i in AA BET BRE CEN ODO OIR SAU JAL RIS ; do mkdir -p $z/$i/IM $z/$i/AM $z/$i/SC $z/$i/IM2M $z/$i/AM2M $z/$i/SC2M $z/$i/SI ; done ; done
#for z in Indep ; do for i in AA BET BRE CEN ODO OIR SAU JAL RIS ; do mkdir -p $z/$i/IM $z/$i/AM $z/$i/SC $z/$i/IM2M $z/$i/AM2M $z/$i/SC2M $z/$i/SI ; done ; done

mv *INDEP.fs* Indep/
mv *_indep.fs* Indep/


mv *Li.fs* Lie/
#Work on INDEP RESULTS:
cd Indep
for j in IM IM2M AM AM2M SC SC2M SI; do mv AA.*_masked_2015_* AA/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv AA_28_28_indep.fs.$j.*_masked_2015_* AA/$j/ ; done 
                                 SI
for j in IM IM2M AM AM2M SC SC2M SI; do mv BET_*INDEP.fs.$j.*_masked_2015_* BET/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv BRE_*INDEP.fs.$j.*_masked_2015_* BRE/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv CEN_*INDEP.fs.$j.*_masked_2015_* CEN/$j/ ; done 
#for j in IM IM2M AM AM2M SC SC2M SI; do mv JAL_*INDEP.fs.$j.*_masked_2015_* JAL/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv ODO_*INDEP.fs.$j.*_masked_2015_* ODO/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv OIR_*INDEP.fs.$j.*_masked_2015_* OIR/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv RIS_*INDEP.fs.$j.*_masked_2015_* RIS/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv SAU_*INDEP.fs.$j.*_masked_2015_* SAU/$j/ ; done 

#Rename folder




cd ../Lie
for j in IM IM2M AM AM2M SC SC2M SI; do mv AA_*LIE.fs.$j.*_masked_2015_* AA/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv BET_*LIE.fs.$j.*_masked_2015_* BET/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv BRE_*LIE.fs.$j.*_masked_2015_* BRE/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv CEN_*LIE.fs.$j.*_masked_2015_* CEN/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv JAL_*LIE.fs.$j.*_masked_2015_* JAL/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv ODO_*LIE.fs.$j.*_masked_2015_* ODO/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv OIR_*LIE.fs.$j.*_masked_2015_* OIR/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv RIS_*LIE.fs.$j.*_masked_2015_* RIS/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv SAU_*LIE.fs.$j.*_masked_2015_* SAU/$j/ ; done 

#Rename folder
cd AA
for i in AM IM IM2M SC SC2M AM2M SI ; do for f in $i/AA_46_46LIE.fs.$i.*_masked_2015_* ; do mv $f ${f/AA_46_46LIE.fs/AA} ; done ; done
cd ../BET
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/BET_24_24LIE.fs.$i.*_masked_2015_* ; do mv $f ${f/BET_24_24LIE.fs/BET} ; done ; done
cd ../BRE
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/BRE_24_24LIE.fs.$i.*_masked_2015_* ; do mv $f ${f/BRE_24_24LIE.fs/BRE} ; done ; done
cd ../CEN
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/CEN_28_28LIE.fs.$i.*_masked_2015_* ; do mv $f ${f/CEN_28_28LIE.fs/CEN} ; done ; done
cd ../JAL
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/JAL_16_16LIE.fs.$i.*_masked_2015_* ; do mv $f ${f/JAL_16_16LIE.fs/CEN} ; done ; done
cd ../ODO
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/ODO_34_34LIE.fs.$i.*_masked_2015_* ; do mv $f ${f/ODO_34_34LIE.fs/ODO} ; done ; done
cd ../OIR/
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/OIR_42_42LIE.fs.$i.*_masked_2015_* ; do mv $f ${f/OIR_42_42LIE.fs/OIR} ; done ; done
cd ../RIS/
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/RIS_28_28LIE.fs.$i.*_masked_2015_* ; do mv $f ${f/RIS_28_28LIE.fs/RIS} ; done ; done
cd ../SAU/
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/SAU_28_28LIE.fs.$i.*_masked_2015_* ; do mv $f ${f/SAU_28_28LIE.fs/SAU} ; done ; done
cd ..

for z in AA BET  OIR RIS  ; do for j in IM IM2M AM AM2M SC SC2M; do  for i in $z/$j/$z* ; do  tail -n 8 $i/*.txt >> $z.$j.txt ; done;  done ; done

#Reshape the results to make it readable in R
for i in AA BET OIR RIS ; do for j in IM IM2M AM AM2M SC SC2M; do grep -v "Grid" $i.$j.txt | grep -E  -A 1 ',|\)]' | grep -v "likelihood" | sed -e 's/--//g' | sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done
#for i in AA BET OIR RIS  ; do for j in IM IM2M AM AM2M SC SC2M ; do grep -v "Grid" $i.$j.txt | grep "," | sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done
for i in AA BET OIR RIS  ; do for j in IM IM2M SC SC2M AM AM2M ; do grep "AIC" $i.$j.txt >> $i.$j.AIC.txt & grep "Optimized log-likelihood" $i.$j.txt >> $i.$j.LogLik.txt & grep "theta" $i.$j.txt >> $i.$j.theta.txt; done ; done
for i in AA BET OIR RIS  ; do for j in IM IM2M SC SC2M AM AM2M ; do paste -d " " $i.$j.AIC.txt $i.$j.LogLik.txt $i.$j.theta.txt $i.$j.Opt.Param.txt >> $i.$j.DADI_1.txt  ; done ; done
for i in AA BET OIR RIS  ; do for j in IM IM2M SC SC2M AM AM2M ; do sed -e 's/AIC://g' $i.$j.DADI_1.txt | sed -e 's/theta://g' | sed -e 's/Optimized log-likelihood://g' | sed -e 's/ptimized parameters: array//g' | sed -e 's/\]//g' | sed -e 's/\[//g' | sed -e 's/(\|)//g' | sed -e '/nan/d' >> $i.$j.DADI_2.txt ; done ; done


sed -e 's/://g' 



#Keep only the "BFGS" results :
# Here: sould use "IF" to treat the txt only if the chain "BFGS" exists in *.txt , otherwise, the file without this chain is siply not considered
#structure: if  "BFGS" exist in *.txt then tail -n... fi 
for z in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M AM AM2M SC SC2M SI; do  for i in $z/$j/$z.$j.* ; do  tail -n 8 $i/*.txt >> $z.$j.txt ; done;  done ; done

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

#####################PEMA#######################################################

for j in IM IM2M AM AM2M SC SC2M SI; do mv AA*PEMA.fs.$j.*_masked_2015_* AA/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv BET_*PEMA.fs.$j.*_masked_2015_* BET/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv BRE_*.fs.$j.*_masked_2015_* BRE/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv CEN_*.fs.$j.*_masked_2015_* CEN/$j/ ; done 
#for j in IM IM2M AM AM2M SC SC2M SI; do mv JAL_*.fs.$j.*_masked_2015_* JAL/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv ODO_*.fs.$j.*_masked_2015_* ODO/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv OIR_*.fs.$j.*_masked_2015_* OIR/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv RIS_*.fs.$j.*_masked_2015_* RIS/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv SAU_*.fs*$j.*_masked_2015_* SAU/$j/ ; done 


for j in IM IM2M AM AM2M SC SC2M SI; do mv  AA.$j.PEMA_masked_2015_* AA/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv BET.$j.PEMA_masked_2015_* BET/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv OIR.$j.PEMA_masked_2015_* OIR/$j/ ; done 
for j in IM IM2M AM AM2M SC SC2M SI; do mv RIS.$j.PEMA_masked_2015_* RIS/$j/ ; done 

#Rename folder
cd AA
for i in AM IM IM2M SC SC2M AM2M ; do for f in $i/AA_24_24_PEMA.fs.$i.*_masked_2015_* ; do mv $f ${f/AA_24_24_PEMA.fs/AA} ; done ; done
cd ../BET                        
for i in AM IM IM2M SC SC2M AM2M ; do for f in $i/BET_24_24PEMA.fs.$i.*_masked_2015_* ; do mv $f ${f/BET_24_24PEMA.fs/BET} ; done ; done
cd ../BRE                        
for i in AM IM IM2M SC SC2M AM2M ; do for f in $i/BRE_24_24PEMA.fs.$i.*_masked_2015_* ; do mv $f ${f/BRE_24_24PEMA.fs/BRE} ; done ; done
cd ../CEN                        
for i in AM IM IM2M SC SC2M AM2M ; do for f in $i/CEN_24_24PEMA.fs.$i.*_masked_2015_* ; do mv $f ${f/CEN_24_24PEMA.fs/CEN} ; done ; done
#cd ../JAL                       
#for i in AM IM IM2M SC SC2M AM2M; do for f in $i/JAL_16_16_PEMA.fs.$i.*_masked_2015_* ; do mv $f ${f/JAL_16_16_PEMA.fs/CEN} ; done ; done
cd ../ODO                        
for i in AM IM IM2M SC SC2M AM2M ; do for f in $i/ODO_34_34PEMA.fs.$i.*_masked_2015_* ; do mv $f ${f/ODO_34_34PEMA.fs/ODO} ; done ; done
cd ../OIR/                       
for i in AM IM IM2M SC SC2M AM2M ; do for f in $i/OIR_30_30PEMA.fs.$i.*_masked_2015_* ; do mv $f ${f/OIR_30_30PEMA.fs/OIR} ; done ; done
cd ../RIS/                       
for i in AM IM IM2M SC SC2M AM2M ; do for f in $i/RIS_20_20PEMA.fs.$i.*_masked_2015_* ; do mv $f ${f/RIS_20_20PEMA.fs/RIS} ; done ; done
cd ../SAU/                       
for i in AM IM IM2M SC SC2M AM2M ; do for f in $i/SAU_28_28PEMA.fs.$i.*_masked_2015_* ; do mv $f ${f/SAU_28_28PEMA.fs/SAU} ; done ; done
cd ..


for z in AA BET BRE CEN ODO OIR RIS SAU ; do for j in SI IM IM2M AM AM2M SC SC2M ; do  for i in $z/$j/$z.$j.* ; do  tail -n 8 $i/*.txt >> $z.$j.txt ; done;  done ; done
#for z in AA BET BRE CEN ODO OIR RIS SAU ; do for j in SI ; do  for i in $z/$j/$z.$j.* ; do  tail -n 8 $i/*.txt >> $z.$j.txt ; done;  done ; done

#Reshape the results to make it readable in R
#for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M AM AM2M SC SC2M; do grep -v "Grid" $i.$j.txt | grep -E  ", | \)]"| sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done
#for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M AM AM2M SC SC2M; do grep -v "Grid" $i.$j.txt | grep -E ',| \)]' | sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done

for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M AM AM2M SC SC2M; do grep -v "Grid" $i.$j.txt | grep -E  -A 1 ',|\)]' | grep -v "likelihood" | sed -e 's/--//g' | sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done

for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M SC SC2M AM AM2M ; do grep "AIC" $i.$j.txt >> $i.$j.AIC.txt & grep "Optimized log-likelihood" $i.$j.txt >> $i.$j.LogLik.txt & grep "theta" $i.$j.txt >> $i.$j.theta.txt; done ; done
for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M SC SC2M AM AM2M ; do paste -d " " $i.$j.AIC.txt $i.$j.LogLik.txt $i.$j.theta.txt $i.$j.Opt.Param.txt >> $i.$j.DADI_1.txt  ; done ; done
for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M SC SC2M AM AM2M ; do sed -e 's/AIC://g' $i.$j.DADI_1.txt | sed -e 's/theta://g' | sed -e 's/Optimized log-likelihood://g' | sed -e 's/ptimized parameters: array//g' | sed -e 's/\]//g' | sed -e 's/\[//g' | sed -e 's/(\|)//g' | sed -e '/nan/d' >> $i.$j.DADI_2.txt ; done ; done
#for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M SC SC2M AM AM2M ; do wc -l $i.$j.DADI_2.txt > $i.$j.Line ; done ; done
#for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M SC SC2M AM AM2M SI; do for `seq in 1 $i.$j.Line` ; do echo $i.$j ; done ; done ; done
rm *Param* *AIC* *DADI_1.txt *theta* *LogLik.txt

mkdir DADI_SYNTHESE_INDEP 
mv *DADI_2.txt DADI_SYNTHESE_INDEP
cd DADI_SYNTHESE_INDEP
#for i in *.DADI_LIE2.txt ; do mv $i ${i%DADI_LIE2.txt}DADI_unfolded.txt ; done

#Compute mean over all files
$libR/dadi_stat_indep_lie.R 

####
paste model DADI.mean.tmp >> mean.tmp.txt
cat entete mean.tmp.txt >> meanFull_INDEP.txt
rm model* entete *tmp*



























libR=/home/ubuntu/Desktop/DADI_RESULT/SCRIPT_ANLAYSE_DADI

#Keep only the "BFGS" results :
# Here: sould use "IF" to treat the txt only if the chain "BFGS" exists in *.txt , otherwise, the file without this chain is siply not considered
#structure: if  "BFGS" exist in *.txt then tail -n... fi 
for z in OIR ; do for j in IM IM2M AM AM2M SC SC2M ; do  for i in /$j/$z.$j_masked* ; do  tail -n 8 $i/*.txt >> $z.$j.txt ; done;  done ; done

for i in OIR ; do for j in IM IM2M AM AM2M SC SC2M; do grep -v "Grid" $i.$j.txt | grep -E  -A 1 ", | \)]"| grep -v "likelihood" | sed -e 's/--//g' | sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done

for i OIR  ; do for j in IM IM2M SC SC2M AM AM2M SI; do grep "AIC" $i.$j.txt >> $i.$j.AIC.txt & grep "Optimized log-likelihood" $i.$j.txt >> $i.$j.LogLik.txt & grep "theta" $i.$j.txt >> $i.$j.theta.txt; done ; done
for i OIR  ; do for j in IM IM2M SC SC2M AM AM2M SI; do paste -d " " $i.$j.AIC.txt $i.$j.LogLik.txt $i.$j.theta.txt $i.$j.Opt.Param.txt >> $i.$j.DADI_1.txt  ; done ; done
for i OIR  ; do for j in IM IM2M SC SC2M AM AM2M SI; do sed -e 's/AIC://g' $i.$j.DADI_1.txt | sed -e 's/theta://g' | sed -e 's/Optimized log-likelihood://g' | sed -e 's/ptimized parameters: array//g' | sed -e 's/\]//g' | sed -e 's/\[//g' | sed -e 's/(\|)//g' | sed -e '/nan/d' >> $i.$j.DADI_2.txt ; done ; done
#for i OIR ; do for j in IM IM2M SC SC2M AM AM2M ; do wc -l $i.$j.DADI_2.txt > $i.$j.Line ; done ; done
#for i OIR ; do for j in IM IM2M SC SC2M AM AM2M SI; do for `seq in 1 $i.$j.Line` ; do echo $i.$j ; done ; done ; done
rm *Param* *AIC* *DADI_1.txt *theta* *LogLik.txt

mkdir DADI_SYNTHESE_INDEP 
mv *DADI_2.txt DADI_SYNTHESE_INDEP
cd DADI_SYNTHESE_INDEP
#for i in *.DADI_LIE2.txt ; do mv $i ${i%DADI_LIE2.txt}DADI_unfolded.txt ; done

#Compute mean over all files
$libR/dadi_stat_indep_lie.R 

####
paste model DADI.mean.tmp >> mean.tmp.txt
cat entete mean.tmp.txt >> meanFull_INDEP.txt
rm model* entete *tmp* DADI.mean.txt








for i in AM  IM2M SC SC2M; do for f in $i/OIR_42_42INDEP.fs.$i.*_masked_2015_* ; do mv $f ${f/OIR_42_42INDEP.fs/OIR} ; done ; done
for i in OIR  ; do for j in IM2M AM SC SC2M ; do grep -v "Grid" $i.$j.txt | grep -E  -A 1 ',|\)]' | grep -v "likelihood" | sed -e 's/--//g' | sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done
for i in OIR  ; do for j in IM2M SC SC2M AM  ; do grep "AIC" $i.$j.txt >> $i.$j.AIC.txt & grep "Optimized log-likelihood" $i.$j.txt >> $i.$j.LogLik.txt & grep "theta" $i.$j.txt >> $i.$j.theta.txt; done ; done
for i in OIR  ; do for j in IM2M SC SC2M AM  ; do paste -d " " $i.$j.AIC.txt $i.$j.LogLik.txt $i.$j.theta.txt $i.$j.Opt.Param.txt >> $i.$j.DADI_1.txt  ; done ; done
for i in OIR  ; do for j in IM2M SC SC2M AM  ; do sed -e 's/AIC://g' $i.$j.DADI_1.txt | sed -e 's/theta://g' | sed -e 's/Optimized log-likelihood://g' | sed -e 's/ptimized parameters: array//g' | sed -e 's/\]//g' | sed -e 's/\[//g' | sed -e 's/(\|)//g' | sed -e '/nan/d' >> $i.$j.DADI_2.txt ; done ; done



cd AA
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/AA_28_28_indep.fs.$i.*_masked_2015_* ; do mv $f ${f/AA_28_28_indep.fs/AA} ; done ; done
cd ../BET                        
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/BET_22_22INDEP.fs.$i.*_masked_2015_* ; do mv $f ${f/BET_22_22INDEP.fs/BET} ; done ; done
cd ../BRE                        
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/BRE_24_24INDEP.fs.$i.*_masked_2015_* ; do mv $f ${f/BRE_24_24INDEP.fs/BRE} ; done ; done
cd ../CEN                        
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/CEN_24_24INDEP.fs.$i.*_masked_2015_* ; do mv $f ${f/CEN_24_24INDEP.fs/CEN} ; done ; done
cd ../JAL                        
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/JAL_16_16INDEP.fs.$i.*_masked_2015_* ; do mv $f ${f/JAL_16_16INDEP.fs/CEN} ; done ; done
cd ../ODO                        
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/ODO_34_34INDEP.fs.$i.*_masked_2015_* ; do mv $f ${f/ODO_34_34INDEP.fs/ODO} ; done ; done
cd ../OIR/                       
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/OIR_30_30INDEP.fs.$i.*_masked_2015_* ; do mv $f ${f/OIR_30_30INDEP.fs/OIR} ; done ; done
cd ../RIS/                       
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/RIS_20_20INDEP.fs.$i.*_masked_2015_* ; do mv $f ${f/RIS_20_20INDEP.fs/RIS} ; done ; done
cd ../SAU/                       
for i in AM IM IM2M SC SC2M AM2M SI; do for f in $i/SAU_28_28INDEP.fs.$i.*_masked_2015_* ; do mv $f ${f/SAU_28_28INDEP.fs/SAU} ; done ; done
cd ..
#Keep only the "BFGS" results :
# Here: sould use "IF" to treat the txt only if the chain "BFGS" exists in *.txt , otherwise, the file without this chain is siply not considered
#structure: if  "BFGS" exist in *.txt then tail -n... fi 
for z in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M AM AM2M SC SC2M SI; do  for i in $z/$j/$z.$j.* ; do  tail -n 8 $i/*.txt >> $z.$j.txt ; done;  done ; done
for z in AA  ; do for j in IM IM2M AM AM2M SC SC2M SI; do  for i in $z/$j/$z.$j_* ; do  tail -n 8 $i/*.txt >> $z.$j.txt ; done;  done ; done
for z in BET BRE CEN ODO OIR RIS SAU  ; do for j in SI; do  for i in $z/$j/$z.$j_* ; do  tail -n 8 $i/*.txt >> $z.$j.txt ; done;  done ; done
for z in OIR  ; do for j in IM2M  AM SC SC2M ; do  for i in $z/$j/$z.$j_* ; do  tail -n 8 $i/*.txt >> $z.$j.txt ; done;  done ; done

#for z in AA BET BRE CEN ODO OIR RIS SAU ; do for j in SI ; do  for i in $z/$j/$z.$j.* ; do  tail -n 8 $i/*.txt >> $z.$j.txt ; done;  done ; done
#for z in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M AM AM2M SC SC2M SI; do  for i in $z/$j/* ; do  tail -n 8 $i >> $z.$j.txt ; done;  done ; done

#Reshape the results to make it readable in R
#for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M AM AM2M SC SC2M; do grep -v "Grid" $i.$j.txt | grep -E  ", | \)]"| sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done
#for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M AM AM2M SC SC2M; do grep -v "Grid" $i.$j.txt | grep -E ',| \)]' | sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done

for i in AA BET BRE CEN ODO OIR RIS SAU ; do for j in IM IM2M AM AM2M SC SC2M SI; do grep -v "Grid" $i.$j.txt | grep -E  -A 1 ',|\)]' | grep -v "likelihood" | sed -e 's/--//g' | sed -e ':z;N;s/\n/ /;bz' | sed -e 's/\(O\|B\)/\n/g' | sed -e 's/\(O\|B\)/\n/g' | sed -e '1d' | sed -e 's/,//g'> $i.$j.Opt.Param.txt ; done ; done

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

####
paste model DADI.mean.tmp >> mean.tmp.txt
cat entete mean.tmp.txt >> meanFull_INDEP.txt
rm model* entete *tmp* DADI.mean.txt


$libR/dadi_stat_indep_lieMIN.R 

####
paste model DADI.min.tmp >> min.tmp.txt
cat entete min.tmp.txt >> minFull_INDEP.txt
rm model* entete *tmp* DADI.mean.txt
