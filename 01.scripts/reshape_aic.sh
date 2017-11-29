tail -n 1 *masked*/*txt |grep -B 1 "AIC" |perl -pe 's/<==\n/\t/g' |less
