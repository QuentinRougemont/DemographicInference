
#very small script to extract AIC of each run after the final BFGS
tail -n 1 *masked*/*txt |\
	grep -B 1 "AIC" |\
	perl -pe 's/<==\n/\t/g' |\
	grep -v "-" | \
	sed 's/:/\t/g' | \
	LC_ALL=C sort -k 4 -n > sorted.models.txt

#list all comparison that have successfully finished:
cut -d "/" -f 1 sorted.models.txt |\
	sort |\
	cut -d "." -f 1 |\
	uniq -c > list_comparaison

#move models:
mv AM* IM* SC* SI* 04-results/
