
for i in AM* IM* SC* SI* ; do
        tail -n 1 $i/*/*txt |\
        grep -B 1 "AIC" |\
        perl -pe 's/<==\n/\t/g' |\
        grep -v "-" | \
        sed 's/:/\t/g' |\
        LC_ALL=C sort -k 4 -n >> models.sorted.tmp ;
done

filename=$( basename $(echo $(pwd)))
sort -k 4,4n -g models.sorted.tmp |\
        sed 's/==>//g' |\
        sed 's/AIC\t//g' > "$filename".sorted.models.txt
rm models.sorted.tmp

ls */*png |cut -d "/" -f 1 |\
        sort |\
        uniq |\
        cut -d "_" -f 1 |\
        sort |\
        uniq -c > "$filename".list_comparaison


