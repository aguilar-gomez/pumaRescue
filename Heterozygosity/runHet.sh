
module load htslib
module load vcftools

#Quick vcftools genome-wide estimate
VCF=bricei_simplePASS_variants_chr.vcf.gz
vcftools --gzvcf $VCF --het --out genome-wide.hetPlink

#Use Jacqueline's vaquita script
tabix -p vcf $VCF
 awk -v OFS="\t" '{print "chr" FNR,$3}' autosomes > chromosome_lengths.txt

#Run script 
#Modify line to make it work on python3:
python3 ./SlidingWindowHet.py bricei_simplePASS_all_$chr.vcf.gz chromosome_lengths.txt 1000000 100000 $chr 