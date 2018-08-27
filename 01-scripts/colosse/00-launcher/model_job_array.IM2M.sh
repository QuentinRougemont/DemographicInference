#!/bin/bash
#PBS -A ihv-653-aa
#PBS -N Output
##PBS -o Out.out
##PBS -e Out.err
#PBS -l walltime=28:00:00
#PBS -l nodes=1:ppn=8
#PBS -M yourmail@yourmail
##PBS -m ea 
#PBS -t [1-4]%4

source /clumeq/bin/enable_cc_cvmfs
source /rap/ihv-653-ab/quentin/01.laben/DemographicInference/temp/bin/activate
# Move to directory where job was submitted
cd "${PBS_O_WORKDIR}"

# Folder to run simulations
fs="your_spectrum.fs" #$1     #dataset fs stored in 03-data
model="IM2m" #$2  #model names
folded=folded #$3 #either "folded" or "unfolded"
grid_size=40 #number of grid points (n pts) in dadi will be n,n+10,n+20
FOLDER=$model.$MOAB_JOBARRAYINDEX

for i in $(seq 8)
do
    ./01-scripts/00.run_dadi_parallel_v2.sh "$fs" "$model" "$folded" "$grid_size" "$FOLDER"_"$i" 
done
