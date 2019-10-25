#!/bin/bash
source /home/quentin/software/DemographicInference/temp/bin/activate
id=$1
fs=$2     #dataset fs stored in 03-data
model=$3  #model name
folded=$4 #either "folded" or "unfolded"
grid_size=$5 #number of grid points (n pts) in dadi will be n,n+10,n+20
grid_size=$(( grid_size   ))
grid_size2=$(( grid_size + 20 ))
grid_size3=$(( grid_size + 30 ))

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
if [ $folded = "folded" ]
then
	python  02-modifs_v2/pang_models.py \
		-o "$model"_"$id" \
		-p $grid_size,$grid_size2,$grid_size3  \
		-f "$fs" \
		-m "$model" \
		-l -v  \
		&>> ./10-log/"$model"_"$id".log;
else
	python  02-modifs_v2/pang_models.py \
	    -o "$model"_$id \
	    -p $grid_size,$grid_size2,$grid_size3  \
	    -f "$fs" \
	    -m "$model" \
	    -l -v \
	    &>> ./10-log/"$model"_"$id".log;
fi 
