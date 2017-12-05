#!/bin/bash

best=$1

name=$(basename $best)

grep -A 7 "BFGS" "$best" |grep -v "BFGS" |sed -e 's/Optimized parameters: array(\[//g' -e 's/Optimized log-likelihood: //g  ' -e 's/theta: //g' -e 's/AIC: //g' -e 's/\])//g' -e 's/,/\n/g' -e 's/ //g' | sed -e '/^$/d' > parameters_best_"${name%.txt}"

