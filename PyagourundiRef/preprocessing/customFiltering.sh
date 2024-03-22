### FilterVCF
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o snpCleaner.log.$JOB_ID.$TASK_ID
#$ -l highp,h_rt=24:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 1
#$ -M daguilar
#$ -t 1-237:1

. /u/local/Modules/default/init/modules.sh

module load python/2.7.15 
module load htslib

IDX=$(printf %03d ${SGE_TASK_ID})
VCF=puma_allsamples_${IDX}_snpEff.vcf.gz 
SCRIPT=~/project-kirk-bigdata/Pconcolor/scripts/filter_puma.py

python2.7 ${SCRIPT} ${VCF} | bgzip > ${VCF%.vcf.gz}_filter.vcf.gz
