### FilterVCF
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o simplify.log.$JOB_ID.$TASK_ID
#$ -l highp,h_rt=24:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 1
#$ -M daguilar

. /u/local/Modules/default/init/modules.sh
module load bcftools
module load htslib

### Simplify & exclude sites failing filters, extract SNPs 
export VCF=puma_simplePASS_SIFT_ALL.vcf.gz

zcat $VCF\
| awk '{ if ( $0~/^#/ || $8~"VariantType=SNP" ){print $0}}' \
| bgzip > puma_simple_SIFT_SNPs.vcf.gz

tabix -p vcf puma_simple_SIFT_SNPs.vcf.gz




