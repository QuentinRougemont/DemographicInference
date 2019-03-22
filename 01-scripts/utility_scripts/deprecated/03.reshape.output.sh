#!/bin/bash
#very basic and not very portable script to reshape output
#you'll probably have to modify it before use

ls -d *masked* > list_file

for i in $(cat list_file) 
do 
    grep -A 9 "BFGS" $i/$i.txt | grep ",\|)"| sed -e ':z;N;s/\n/ /;bz' -e 's/,//g' |sed -e 's/Optimized parameters: array//g'  -e 's/\]//g' -e 's/\[//g' -e 's/,//g' -e 's/(//g' -e 's/)//g' > 04-results/$i.Opt.Param.txt ; 
    grep -A 9 "BFGS" $i/$i.txt | grep "AIC" | sed -e 's/AIC: //g' > 04-results/$i.AIC.txt  
    grep -A 9 "BFGS" $i/$i.txt | grep "Optimized log-likelihood" |sed -e 's/Optimized log-likelihood: //g'> 04-results/$i.LogLik.txt
    grep -A 9 "BFGS" $i/$i.txt | grep "theta" | sed -e 's/theta: //g' > 04-results/$i.theta.txt; 
done

for i in $(cat list_file) 
do 
    paste -d " " 04-results/$i.AIC.txt 04-results/$i.LogLik.txt 04-results/$i.theta.txt 04-results/$i.Opt.Param.txt | sed -e "s/^/$i\t/g"  > 04-results/$i.DADI_1.txt 
done

cat 04-results/SC_*.DADI_1.txt |LC_ALL=C sort -k 2 > SC.DADI.txt
cat 04-results/SC2M_*.DADI_1.txt |LC_ALL=C sort -k 2  >SC2M.DADI.txt
cat 04-results/IM_*.DADI_1.txt |LC_ALL=C sort -k 2  > IM.DADI.txt
cat 04-results/IM2M_*.DADI_1.txt |LC_ALL=C sort -k 2  > IM2M.DADI.txt
cat 04-results/AM_*.DADI_1.txt |LC_ALL=C sort -k 2 > AM.DADI.txt
cat 04-results/AM2M_*.DADI_1.txt |LC_ALL=C sort -k 2  >AM2M.DADI.txt
cat 04-results/SI_*.DADI_1.txt |LC_ALL=C sort -k 2  > SI.DADI.txt

cat *.DADI.txt |LC_ALL=C sort -k 2 > ALL.models.txt 
