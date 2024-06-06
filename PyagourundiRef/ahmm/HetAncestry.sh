#Concatenate all files of heterozygosity
cat *txt|grep -v "chrom" > allGenome.Het_noheader
head -1 puma_NW_024412376.1.vcf.gz_het_1000kbWin_100Kbstep.txt > header
#Add header back
cat header allGenome.Het_noheader |cut -f1-2,54- > AllHeterozygosity.txt
#Add column to convert to bedfile format
awk 'BEGIN {OFS="\t"} {$3 = $2 + 1000000;$2 = $2 - 1;  print}' AllHeterozygosity.txt > Heterozygosity.bed

#Convert ahmm outout into bedfile format
awk 'BEGIN {OFS="\t"} {$4 =int($3 * 2) / 2;$3 = $2;$2 = $2 - 1;print}' AFP8.yag.prob.filter0.95|tail -n+2  >AFP8.yag.prob.filter0.95.bed
#Extract colum with AFP8 heterozygosity
cut -f1-3,32 Heterozygosity.bed|tail -n+2 > AFP8.het.bed
#bedtools intersect to get windows with ancestry and heterozygosity
bedtools intersect -a AFP8.yag.prob.filter0.95.bed -b AFP8.het.bed  -wb > AFP8.HetandAncestry.bed



	
