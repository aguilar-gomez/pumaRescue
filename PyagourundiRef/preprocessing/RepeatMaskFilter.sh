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


### VariantFiltration to mask repeats
# Note: GATK requires indexed mask file
gatk IndexFeatureFile -I ${REPEATMASK}
#gatk IndexFeatureFile 

IDX=$(printf %03d ${SGE_TASK_ID})
VCF=bricei_joint_${IDX}_snpEff_SIFT.vcf.gz  
REGION=$(ls $(dirname ${REFERENCE})/intervals/*_${IDX}.bed)

gatk VariantFiltration \
-R ${REFERENCE} \
-L ${REGION} \
-mask ${REPEATMASK} --mask-name "FAIL_Repeat" \
-verbosity ERROR \
-V $VCF \
-O ${VCF%.vcf.gz}_LeftAlignTrim_Mask.vcf.gz


### Custom filtering (see bricei_filter.py) 

#Get mean coverge from DepthOfCoverage run
cut -f1,3 -d "," *sample_summary| grep "Bede" > meanCovBricei

#Filter

module load htslib

export SCRIPT=~/project-klohmuel/Bricei/scripts/bricei_filter.py
IDX=$(printf %03d ${SGE_TASK_ID})
export VCF=bricei_joint_${IDX}_snpEff_SIFT_LeftAlignTrim_Mask.vcf.gz  

python2.7 ${SCRIPT} ${VCF} | bgzip > ${VCF%.vcf.gz}_Filter.vcf.gz

tabix -p vcf ${VCF%.vcf.gz}_Filter.vcf.gz


### Simplify & exclude sites failing filters, extract SNPs 
IDX=$(printf %03d ${SGE_TASK_ID})
export VCF=bricei_joint_${IDX}_snpEff_SIFT_LeftAlignTrim_Mask_Fil.vcf.gz

bcftools annotate \
-x ^INFO/AC,INFO/AF,INFO/AN,INFO/ANN,INFO/LOF,INFO/NMD,INFO/SIFTINFO,INFO/VariantType,FORMAT \
-Ov $VCF \
| awk '{ if ( $0~/^#/ || $7=="PASS" ){print $0}}' \
| bgzip > bricei_${IDX}_simple_PASS.vcf.gz

tabix -p vcf bricei_${IDX}_simple_PASS.vcf.gz

zcat bricei_${IDX}_simple_PASS.vcf.gz\
| awk '{ if ( $0~/^#/ || $8~"VariantType=SNP" ){print $0}}' \
| bgzip > bricei_${IDX}_simple_PASS_variants.vcf.gz

tabix -p vcf bricei_${IDX}_simple_PASS_variants.vcf.gz

### Concatenate

ls -v bricei_*_simple_PASS.vcf.gz > simplePass.vcflist
ls -v bricei_*_simple_PASS_variants.vcf.gz > Variants.vcflist

export NUMTHREADS=8

bcftools concat -f simplePass.vcflist --threads ${NUMTHREADS} -Oz -o bricei_simplePASS_all.vcf.gz
bcftools concat -f variants.vcflist --threads ${NUMTHREADS} -Oz -o bricei_simplePASS_variants_all.vcf.gz


#Keep only autosomes
grep "CM" *fai > chromosomes.fai
grep -v CM051429.1 chromosomes.fai > autosomes.fai
awk -v OFS="\t" '{print $1, 1, $2}' autosomes.fai > autosomes

bcftools index bricei_simplePASS_variants_all.vcf.gz
bcftools view -R autosomes bricei_simplePASS_variants_all.vcf.gz -Oz -o bricei_simplePASS_variants_autosomes.vcf.gz

#Rename chrs
awk -v OFS="\t" '{print $1, "chr" FNR}' autosomes > chromosome_mapping
IN=bricei_simplePASS_variants_autosomes.vcf.gz
OUT=bricei_simplePASS_variants_chr.vcf.gz
bcftools annotate --rename-chrs chromosome_mapping --output $OUT --output-type z $IN
