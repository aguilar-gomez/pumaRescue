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

IDX=$(printf %03d ${SGE_TASK_ID})
VCF=puma_allsamples_${IDX}_snpEff.vcf.gz 

SNPCLEANER=~/project-kirk-bigdata/Pconcolor/scripts/snpCleaner.pl
# SNPcleaner v2.4.3 #
#-v to keep invariant sites too
#-H FLOAT min p-value for exact test of excess of heterozygous [0]
#-h FLOAT min p-value for exact test of HWE [0.0001]
#-b FLOAT min p-value for base quality bias [1e-100]
#-S FLOAT min p-value for strand bias [0.0001]
#-f FLOAT min p-value for map quality bias [0]
#-e FLOAT min p-value for end distance bias [0.0001]
#-d INT   minimum raw site read depth (VCF DP flag) [1]
#-D INT   maximum raw site read depth (VCF DP flag)[1000000]

#mean depth of an individual 19.1 * 50 = 954
#median depth of an individual 12.7 * 50 = 634

#median depth of an individual 12.7 * 50 = 634 * (1/3) = 211
#median depth of an individual 12.7 * 50 = 634 * (2) = 1267

#ultimately calibrated with a histogram of depth in the largest scaffold

$SNPCLEANER  -v -H 1e-6 -h 1e-4 -b 1e-100 \
             -S 1e-4 -f 1e-6 -e 1e-4 -B filtered_${IDX}_all2yag.bed -p all2yag_${IDX}_failedsites.txt.bz \ 
             -d 200 -D 1600 $VCF | bgzip > puma_${IDX}_filter.vcf.gz 
