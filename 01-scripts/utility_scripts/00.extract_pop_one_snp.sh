#!/bin/bash

pop1=$1
pop2=$2
vcffile=$3

if [[ -z "$pop1" ]]
then
       echo "Error: need pop name"
       echo "please provide populations name and vcffile:"
       echo "./00.extract_pop.sh pop1 pop2 yourvcf.vcf"
       exit
fi

if [[ -z "$pop2" ]]
then
       echo "Error: need pop name"
       echo "please provide populations name and vcffile:"
       echo "./00.extract_pop.sh pop1 pop2 yourvcf.vcf"
       exit
fi

if [[ -z "$vcffile" ]]
then
       echo "Error: need input vcf-file!!!"
       echo "please provide populations name and vcffile:"
       echo "./00.extract_pop.sh pop1 pop2 yourvcf.vcf"
       exit
fi

folder="data_"$pop1"_"$pop2""
mkdir "$folder"

#load vcftools if necessary
#module load vcftools 
vcftools --vcf "$vcffile" --keep "$pop1" --max-missing 0.95 --hwe 0.01 --out "$folder"/batch_"$pop1" --recode
vcftools --vcf "$vcffile" --keep "$pop2" --max-missing 0.95 --hwe 0.01 --out "$folder"/batch_"$pop2" --recode

# Compress & index
bgzip "$folder"/batch_"$pop1".recode.vcf
tabix -p vcf "$folder"/batch_"$pop1".recode.vcf.gz
bgzip "$folder"/batch_"$pop2".recode.vcf
tabix -p vcf "$folder"/batch_"$pop2".recode.vcf.gz

# Merge VCF 
vcf-merge "$folder"/batch_"$pop1".recode.vcf.gz "$folder"/batch_"$pop2".recode.vcf.gz > "$folder"/batch_"$pop1"."$pop2".vcf

# Keeps SNPs in the 2 pops
grep -v "^#" "$folder"/batch_"$pop1"."$pop2".vcf | grep -e "0,1" > "$folder"/batch_"$pop1"."$pop2".spectrum.tmp 
grep "^#" "$folder"/batch_"$pop1"."$pop2".vcf > header
cat header "$folder"/batch_"$pop1"."$pop2".spectrum.tmp > "$folder"/batch_"$pop1"."$pop2".spectrum.vcf 
rm header "$folder"/*tmp 

#keep one single SNP by loci
less "$folder"/batch_"$pop1"."$pop2".spectrum.vcf |grep -v "#" |awk '{print $3}' |sed 's/_/\t/g' | awk '{print $2, $1}' |uniq -f 1  |awk '{print $2"_"$1}' > "$folder"/markers
cd "$folder"/
Rscript ../01-scripts/utility_scripts/get_single_snp.R batch_"$pop1"."$pop2".spectrum.vcf markers
cd ../
grep "#" "$folder"/batch_"$pop1"."$pop2".spectrum.vcf > "$folder"/header.tmp

cat "$folder"/header.tmp "$folder"/vcf_one_snp.tmp > "$folder"/batch_"$pop1"."$pop2".spectrum.one.snp.vcf 
rm "$folder"/*tmp

#remove globally invariant sites
vcftools --vcf "$folder"/batch_"$pop1"."$pop2".spectrum.one.snp.vcf --freq --out "$folder"/freq_"$pop1"."$pop2" 
less "$folder"/freq_"$pop1"."$pop2".frq |grep -v ":1" |awk '{print $1, $2}' |sed 's/ /\t/g' > "$folder"/polymorphic_sites 
vcftools --vcf  "$folder"/batch_"$pop1"."$pop2".spectrum.vcf --positions "$folder"/polymorphic_sites --out "$folder"/batch_"$pop1"."$pop2".spectrum.polym  --recode --recode-INFO-all 

#ici on peut rouler le script pour ne garder qu'un seul SNPs par RAD loci
#bash get_single.snp.sh "$folder"/batch_"$pop1"."$pop2".spectrum.vcf 
#Extract again the pop
vcftools --vcf "$folder"/batch_"$pop1"."$pop2".spectrum.polym.recode.vcf --keep "$pop1" --out "$folder"/batch_"$pop1"_spectrum --recode
vcftools --vcf "$folder"/batch_"$pop1"."$pop2".spectrum.polym.recode.vcf --keep "$pop2" --out "$folder"/batch_"$pop2"_spectrum --recode

#allele frequencies
vcftools --vcf "$folder"/batch_"$pop1"_spectrum.recode.vcf --freq --out "$folder"/batch_"$pop1"_spectrum
vcftools --vcf "$folder"/batch_"$pop2"_spectrum.recode.vcf --freq --out "$folder"/batch_"$pop2"_spectrum

#Merge the files
paste "$folder"/batch_"$pop1"_spectrum.frq  "$folder"/batch_"$pop2"_spectrum.frq | perl -pe 's/:/\t/g' | sed '1d' > "$folder"/batch_spectrum.txt
awk '{print $5, $5, $4, $5, $6, $7, $8, $12,  $14, $16, $1, $2 }' "$folder"/batch_spectrum.txt > "$folder"/spectrum.tmp

cd "$folder"/
Rscript ../01-scripts/utility_scripts/01.prepare.spectrum.R 
#To do: replace Rscript by a single bash script
exit
