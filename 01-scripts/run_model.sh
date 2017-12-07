#!/bin/bash

fs=$1     #dataset fs stored in 03-data
model=$2  #model names
folded=$3 #either "folded" or "unfolded"

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

for i in $(seq 1 ); 
do
    if [ $folded = "folded" ]
    then
        python ../02-modifs/script_inference_anneal2_newton_mis_new_models_folded2.py -o "$model"_$i -y pop1 -x pop2 -p 40,50,60 -f "$fs" -m "$model" -l -v  -z &>> ../10-log/"$model"_"$i".log;
    else
        python 02-modifs/unfolded/script_inference_anneal2_newton_mis.py -o 04-results/"$model"_$i -y pop1 -x pop2 -p 40,50,60
        -f "$fs" -m "$model" -l -v  -z &>> 10-log/"$model"_"$i".log;
   fi 
done 

mv *masked* 04-results
