#!/bin/bash

id=$1
fs=$2     #dataset fs stored in 03-data
model=$3  #model names
folded=$4 #either "folded" or "unfolded"
grid_size=$5

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

if [ $folded = "folded" ]
then
	python ../02-modifs_v2/folded/script_inference_demo_new_models_folded.py \
		-o "$model"_"$id" \
		-y pop1 -x pop2 \
		-p $grid_size,$grid_size2,$grid_size3  \
		-f "$fs" \
		-m "$model" \
		-l -v  \
		-z \
		&>> ../10-log/"$model"_"$id".log;
else
    python ../02-modifs_v2/unfolded/script_inference_demo_new_models.py \
	    -o "$model"_$id \
	    -y pop1 -x pop2 \
	    -p $grid_size,$grid_size2,$grid_size3  \
	    -f "$fs" \
	    -m "$model" \
	    -l -v \
	    -z \
	    &>> ./10-log/"$model"_"$id".log;
fi 
