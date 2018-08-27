#!/bin/bash
#source temp/bin/activate
# Move to directory where job was submitted

#List of wanted arguments
fs=$1     #JSFS stored in 03-data
model=$2  #model name
folded=$3 #either "folded" or "unfolded"
grid_size=$4 #number of grid points (n pts) in dadi will be n,n+10,n+20

#do not forget to update the crwd (see dadi manual)
if [[ -z "$fs" ]]
then
    echo "Error: need input file"
    exit
fi

if [[ -z "$model" ]]
then
    echo "Error: need model name ex: SC2M"
    exit
fi

if [[ -z "$folded" ]]
then
    echo "Error: need folding info"
    exit
fi

if [[ -z "$grid_size" ]]
then
    echo "Error: need gride_size for optimisation"
    exit
fi

nrep=32
NUM_CPUS=32
seq $nrep |parallel -j "$NUM_CPUS" ./01-scripts/01-run_model_iteration_v2.sh {} "$fs" "$model" "$folded" "$grid_size"
