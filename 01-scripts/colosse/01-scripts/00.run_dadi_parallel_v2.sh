#!/bin/bash
# Move to directory where job was submitted
#List of wanted arguments
fs=$1     #dataset fs stored in 03-data
model=$2  #model names
folded=$3 #either "folded" or "unfolded"
grid_size=$4
FOLDER=$5

mkdir "$FOLDER" 2>/dev/null
cd "$FOLDER"
cp ../03-data/"$fs" .

nrep=1
NUM_CPUS=1
seq $nrep |parallel -j "$NUM_CPUS" ../01-scripts/01-run_model_iteration_v2.sh {} "$fs" "$model" "$folded" "$grid_size"
