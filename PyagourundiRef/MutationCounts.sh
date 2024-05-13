bcftools view bricei_minGQ30_SNPs_chr.vcf.gz |grep "|SYNONYMOUS|\|#" > bricei_SYNONYMOUS.vcf
bcftools view bricei_minGQ30_SNPs_chr.vcf.gz |grep "|NONSYNONYMOUS|\|#" > bricei_NONSYNONYMOUS.vcf


bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -H bricei_SYNONYMOUS.vcf  > GT_syn
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -H bricei_NONSYNONYMOUS.vcf  > GT_non

/u/home/d/daguilar/project-klohmuel/Bricei/analysis/EmpiricalLoad


#PUMAS

#!/bin/bash
#$ -cwd
#$ -j y
#$ -o concatALL
#$ -l highp,h_rt=24:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 1
#$ -M daguilar


. /u/local/Modules/default/init/modules.sh
module load bcftools
module load htslib

ls -v puma_*_simple_PASS_SIFT.vcf.gz  > list

export NUMTHREADS=8
bcftools concat -f list --threads ${NUMTHREADS} -Oz -o puma_simplePASS_SIFT_ALL.vcf.gz

####################################################################################################

bcftools view puma_simplePASS_SIFT.vcf.gz |grep "|SYNONYMOUS|\|#" > puma_SYNONYMOUS.vcf
bcftools view puma_simplePASS_SIFT.vcf.gz |grep "|NONSYNONYMOUS|\|#" > puma_NONSYNONYMOUS.vcf


#Extract the sites directly
cat SIFT*/*xls|grep "DELETERIOUS" > deleterious.SNPs
cat SIFT*/*xls|grep "TOLERATED" > tolerated.SNPs
cut -f9 tolerated.SNPs |sort |uniq -c
#SIFT results
Deleterious variants
 218,226 NONSYNONYMOUS
    518 START-LOST
   1,680 SYNONYMOUS

Tolerated variants
256,972 NONSYNONYMOUS
    132 START-LOST
 521,403 SYNONYMOUS


bcftools view puma_simplePASS_SIFT_ALL.vcf.gz -R tolerated.SNPs > puma_tolerated.vcf
bcftools view puma_simplePASS_SIFT_ALL.vcf.gz -R deleterious.SNPs > puma_deleterious.vcf

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -H puma_tolerated.vcf  > GT_puma_tol
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -H puma_deleterious.vcf  > GT_puma_del
sed -i '1 s/\[[0-9]*\]//g; 1 s/# //; 1 s/:GT//g' GT_puma_del
sed -i '1 s/\[[0-9]*\]//g; 1 s/# //; 1 s/:GT//g' GT_puma_tol



bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -H puma_SYNONYMOUS.vcf > GT_puma_syn
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -H puma_NONSYNONYMOUS.vcf  > GT_puma_non
sed -i '1 s/\[[0-9]*\]//g; 1 s/# //; 1 s/:GT//g' GT_puma_syn
sed -i '1 s/\[[0-9]*\]//g; 1 s/# //; 1 s/:GT//g' GT_puma_non


bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -H puma_simplePASS_SIFT_ALL.vcf.gz > GT_all_pumas
sed -i '1 s/\[[0-9]*\]//g; 1 s/# //; 1 s/:GT//g' GT_all_pumas

stats --samples '-' puma_simplePASS_SIFT_ALL.vcf.gz --threads 10 |grep "PSC" > counts_PUMAS
sed  '1 s/\[[0-9]*\]//g; 1 s/# //;' counts_PUMAS

grep -v "Note" counts_PUMAS| cut -f 3-6,14 > counts2normalize


#Fixed differences between reference (PYag) and Pconcolor
 bcftools view puma_deleterious.vcf|grep -v "0/"|grep -v "#"|wc
#  2241  132219 3443615
bcftools view puma_tolerated.vcf|grep -v "0/"|grep -v "#"|wc
#  41145 2427555 62501389
