#!/bin/bash

id=$1
fs=$2     #dataset fs stored in 03-data
model=$3  #model name
folded=$4 #either "folded" or "unfolded"
grid_size=$5 #number of grid points (n pts) in dadi will be n,n+10,n+20
grid_size=$(( grid_size * 2 ))
grid_size2=$(( grid_size + 10 ))
grid_size3=$(( grid_size + 20 ))

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

if [[ -z "$grid_size" ]]
then
    echo "Error: need gredd size"
    exit
fi
#run models withou masking
if [[ $folded = "folded" ]]
then
	echo "running folded model"
	echo "by default singleton are not masked (used -z option)"
	python ./02-modifs_v2/folded/script_inference_demo_new_models_folded.py \
		-o "$model"_"$id" \
		-y pop1 -x pop2 \
		-p $grid_size,$grid_size2,$grid_size3  \
		-f "$fs" \
		-m "$model" \
		-l -v  \
		&>> ./10-log/"$model"_"$id".log;
else
	echo "running folded model"
	echo "by default singleton are not masked (used -z option)"
    python ./02-modifs_v2/unfolded/script_inference_demo_new_models.py \
	    -o "$model"_$id \
	    -y pop1 -x pop2 \
	    -p $grid_size,$grid_size2,$grid_size3  \
	    -f "$fs" \
	    -m "$model" \
	    -l -v \
	    &>> ./10-log/"$model"_"$id".log;
fi 
