#!/bin/bash
#PBS -A ihv-653-aa
#PBS -N OutputTest
##PBS -o OutTest.out
##PBS -e OutTest.err
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=8
#PBS -M quentinrougemont@orange.fr
##PBS -t [1-4]%4

source /clumeq/bin/enable_cc_cvmfs
source /rap/ihv-653-ab/quentin/01.laben/DemographicInference/temp/bin/activate
# Move to directory where job was submitted
cd "${PBS_O_WORKDIR}"

# Folder to run simulations
fs="male_120.fs" #$1     #dataset fs stored in 03-data
model="PANGb2N"     #model names
folded=folded #$3 #either "folded" or "unfolded"
grid_size=120 #$4
FOLDER=$model.$MOAB_JOBARRAYINDEX

for i in $(seq 4)
do
    ./01-scripts/00.pan.sh "$fs" "$model" "$folded" "$grid_size" "$FOLDER"_"$i" 
    #./"$MODEL" "$FOLDER"_"$i"  &
done
