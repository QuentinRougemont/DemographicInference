#very small script to extract AIC of each run after the final BFGS
for i in AM* IM* SC* SI* ; do 
	tail -n 1 $i/*masked*/*txt |\
	grep -B 1 "AIC" |\
	perl -pe 's/<==\n/\t/g' |\
	grep -v "-" | \
	sed 's/:/\t/g' |\
	LC_ALL=C sort -k 4 -n >> models.sorted.tmp ; 
done                              

sort -k 4 -n models.sorted.tmp > sorted.models.txt
rm models.sorted.tmp
#list all comparison that have successfully finished:
cut -d "/" -f 1 sorted.models.txt |\
    sort |\
    cut -d "." -f 1 |\
    uniq -c > list_comparaison
#move models:
mv AM* IM* SC* SI* 04-results/
