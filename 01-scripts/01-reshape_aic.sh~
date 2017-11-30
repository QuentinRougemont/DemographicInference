tail -n 1 *masked*/*txt |grep -B 1 "AIC" |perl -pe 's/<==\n/\t/g' |grep -v "-" | sed 's/:/\t/g' |LC_ALL=C sort -k 3 > models.sorted.txt
