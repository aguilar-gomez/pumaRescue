### FilterVCF
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o simplify.log.$JOB_ID.$TASK_ID
#$ -l highp,h_rt=24:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 1
#$ -M daguilar
#$ -t 1-237:1

### Simplify & exclude sites failing filters, extract SNPs 
IDX=$(printf %03d ${SGE_TASK_ID})
export VCF=puma_allsamples_${IDX}_snpEff_filter_LeftAlignTrim_Mask_noSex.vcf.gz 

bcftools annotate \
-x ^INFO/AC,INFO/AF,INFO/AN,INFO/ANN,INFO/LOF,INFO/NMD,INFO/SIFTINFO,INFO/VariantType,FORMAT \
-Ov $VCF \
| awk '{ if ( $0~/^#/ || $7=="PASS" ){print $0}}' \
| bgzip > puma_${IDX}_simple_PASS.vcf.gz

tabix -p vcf puma_${IDX}_simple_PASS.vcf.gz

zcat puma_${IDX}_simple_PASS.vcf.gz\
| awk '{ if ( $0~/^#/ || $8~"VariantType=SNP" ){print $0}}' \
| bgzip > puma_${IDX}_simple_PASS_variants.vcf.gz

tabix -p vcf puma_${IDX}_simple_PASS_variants.vcf.gz

### Concatenate ###################################################################################################

ls -v puma_*_simple_PASS.vcf.gz > simplePass.vcflist
ls -v puma_*_simple_PASS_variants.vcf.gz > Variants.vcflist

export NUMTHREADS=8

bcftools concat -f simplePass.vcflist --threads ${NUMTHREADS} -Oz -o puma_simplePASS_all.vcf.gz
bcftools concat -f variants.vcflist --threads ${NUMTHREADS} -Oz -o puma_simplePASS_variants_all.vcf.gz

