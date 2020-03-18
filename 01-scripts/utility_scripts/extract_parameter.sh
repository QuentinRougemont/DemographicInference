#!/bin/bash

#script to reshape dadi estimates of parameter into one single line
#takes as input the name of the demographic model tested

model=$1  #full name of the model folder
tail -n 6 $model/$model.txt |perl -pe 's/\n/ /g' |\
      sed 's/,/ /g' |\
      sed 's/Optimized parameters: array(\[//g' |\
      sed 's/\]) Optimized log-likelihood://g' |\
      sed "s/ \+/ /g"  |\
      sed 's/[a-zA-Z:]*//g' |\
      awk -v AR="$model" '{ printf "%s",AR OFS  $(NF) OFS $(NF-2) OFS $(NF-1) OFS ; for(i=1; i<NF-2; i++) printf "%s",$i OFS; printf ORS}'
      awk -v AR="$model" '{ printf "%s",AR OFS  $(NF) OFS $(NF-2) OFS $(NF-1) OFS ; for(i=1; i<NF-2; i++) printf "%s",$i OFS; printf ORS}' > param_"$model"
