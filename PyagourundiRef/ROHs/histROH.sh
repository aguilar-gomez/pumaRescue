#Extract individuals for histogram 


#PLINK ROHS
bricei_autosomes_ROH.hom

cut -d " " -f1 *fam > samples

while IFS= read -r sample
do
    number=$(echo "$sample" | tr -dc '0-9')
    echo "$sample $number"
    grep "$sample" bricei_autosomes_ROH.hom | cut -f1,9 > "bricei${number}_rohsKB"
done < "samples"


#VCFTOOLS ROHS
while IFS= read -r sample
do
    number=$(echo "$sample" | tr -dc '0-9')
    echo "$sample $number"
    grep -h "$sample" bricei_simplePASS_variants_chr.vcf.gz_vcftools_chr*.LROH | cut -f1-3 | awk '{print $1 "\t" ($3-$2)/1000}'> "bricei${number}_rohsKB_vcftools"
done < "samples"

#BCFTOOLS ROHS
awk 'BEGIN{OFS="\t"} {print $2, $3, $4, $5, $6, $7, $8}' bricei_simplePASS_variants_autosomes.vcf.gz_bcftoolsROH.txt > ROHs_autosomes_bcftools.tab
while IFS= read -r sample
do
    number=$(echo "$sample" | tr -dc '0-9')
    echo "$sample $number"
    grep "$sample" ROHs_autosomes_bcftools.tab | cut -f1,5 > "bricei${number}_rohs_bcftools"
done < "samples"

