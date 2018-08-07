#!/bin/bash
#source temp/bin/activate
# Move to directory where job was submitted
cd 01-scripts

#List of wanted arguments
fs=$1     #dataset fs stored in 03-data
model=$2  #model names
folded=$3 #either "folded" or "unfolded"
grid_size=$4

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

nrep=64
NUM_CPUS=32
seq $nrep |parallel -j "$NUM_CPUS" ./01-run_model_iteration_v2.sh {} "$fs" "$model" "$folded" "$grid_size"
