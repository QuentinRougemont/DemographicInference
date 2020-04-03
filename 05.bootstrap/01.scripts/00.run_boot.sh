#!/bin/bash

#required arguments: name of the best model
#easy script to construct bootstrap jSFS dataset with ms
#based on dadi parameter estimate
#if on clusters:
#~source /home/qurou/software/DemographicInference/temp/bin/activate

if [ $# -ne 3 ]; then
    echo "USAGE: $0 model P Q"
    echo "Expecting the following values on the command line, in that order"
    echo "Name of the best model to bootstrap (SC2N2mG  or IM2N2mG for now)"
    echo "value of the P parameter"
    echo "value of the P parameter"
    exit 1
else
    model=$1   #name of best model (either SC2N2mG or IM2N2mG  for now)
    P=$2       #value of P parameter in dadi
    Q=$3       #value of Q parameter in dadi
    echo "model is $model"
    echo "P value is $P"
    echo "Q value is $Q"
fi


mkdir simuls_ms_${model}/

#create frequency spectrum for:
#1 "neutral" region
#2 "genomic islands; me'"
#3 "low recbom (linked sel); Ne scaled by Hrf"
#4 "normall recomb"

nboot=100 #create 100 bootstrap dataset (we may create more)
for i in $(seq $nboot) ; do
    01.scripts/boot/boot.ms.${model}_island.sh ${model} >simuls_ms_${model}/data.$i.islands
    01.scripts/boot/boot.ms.${model}_linked_sel.sh ${model} >simuls_ms_${model}/data.$i.linkedsel
    01.scripts/boot/boot.ms.${model}_normal_recomb.sh ${model} >simuls_ms_${model}/data.$i.normal_recomb
    01.scripts/boot/boot.ms.${model}.sh ${model} >simuls_ms_${model}/data.$i.neutral ;
done

mkdir 02.boot_fs_$model
for i in $(seq $nboot) ; do
    ./01.scripts/dadi_create_fs.py simuls_ms_${model}/data.$i. ${P} ${Q} $i ${model} ;
done

