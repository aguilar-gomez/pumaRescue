
git clone --recursive https://github.com/rvaser/sift4g.git sift4
git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git scripts_to_build_SIFT_db
cd test_files/
cp homo_sapiens-test.txt puma_concolor.txt





#!/bin/bash
#$ -cwd
#$ -j y
#$ -o SIFT.log.$JOB_ID.$TASK_ID
#$ -l highp,h_rt=72:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 2
#$ -M daguilar
#$ -t 1-237:1

. /u/local/Modules/default/init/modules.sh


module load java/jdk-11.0.14
module load htslib
module load bcftools


IDX=$(printf %03d ${SGE_TASK_ID})
export SIFT=~/project-klohmuel/programs/SIFT4G_Annotator.jar
export DATABASE=~/project-kirk-bigdata/Pconcolor/genome_outgroup/SIFTDatabse
export VCF=puma_${IDX}_simple_PASS.vcf.gz.tbi
export LOG=${VCF%.vcf.gz}_SIFT.log


echo -e "[$(date "+%Y-%m-%d %T")] Creating working directory..."
DATE=$(date +%Y%m%d%H%M)
TEMPDIR="SIFT${SGE_TASK_ID}_${DATE}"
mkdir $TEMPDIR

gzip -cd ${VCF} > ${VCF%.gz}

java -jar ${SIFT} -c -t -d ${DATABASE} -i ${VCF%.gz} -r  $TEMPDIR

cd $TEMPDIR

#Change format from SIFT to GATK
#the header lines must be formatted with "=" and not ":" and you are not supposed to have spaces within the INFO column
awk '{if ($1~"^#"){ gsub("##SIFT_Threshold: 0.05", "##SIFT_Threshold=0.05"); print } \
else { gsub(" ", "_"); print }}' ${VCF%.vcf.gz}_SIFTpredictions.vcf \
| bgzip > ${VCF%.vcf.gz}_SIFT.vcf.gz

tabix -p vcf ${VCF%.vcf.gz}_SIFT.vcf.gz

rm ../${VCF%.gz}
rm ${VCF%.vcf.gz}_SIFTpredictions.vcf

