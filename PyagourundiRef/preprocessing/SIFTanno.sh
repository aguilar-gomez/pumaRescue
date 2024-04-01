 ### SIFT annotation
# Note: Following instructions at https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB,
# Jacqueline created the database

cp /u/project/klohmuel/DataRepository/Rice_whale/sift_mBalRic1.hap2.tar.gz ~/project-klohmuel/programs/SIFTdatabases/
tar -xzvf ~/project-klohmuel/programs/SIFTdatabases/*.tar.gz
# SIFT4G algorithm (https://github.com/rvaser/sift4g).
### Mask repeats (from RepeatMasker and TandemRepeatsFinder)

module load java/jdk-11.0.14
module load bcftools

export SIFT=~/project-klohmuel/programs/SIFT4G_Annotator.jar
export DATABASE= ~/project-klohmuel/programs/SIFTdatabases/mBalRic1.hap2
export VCF=bricei_joint_${IDX}_snpEff.vcf.gz  
export LOG=${VCF%.vcf.gz}_SIFT.log

echo -e "[$(date "+%Y-%m-%d %T")] Creating working directory..."
DATE=\$(date +%Y%m%d%H%M)
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






