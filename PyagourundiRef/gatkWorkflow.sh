### Bricei calling haplotypes and filtering pipeline (example commands)

# Workflow
# - Index genome and generate chunks
# - Per sample:
#   - Rename bams to Bricei (previously Bede)
#   - DepthOfCoverage
#   - HaplotypeCaller
# - All samples jointly:
#   - GenotypeGVCFs 
#   - LeftAlignAndTrimVariants
#   - SnpEff
#   - SIFT
#   - Mask repeats (identified with RepeatMasker + Tandem Repeats Finder) 
#   - Custom filtering
#   - Simplify, extract SNPs
#   - Concatenate VCF files
#   - Extract autosomes
################################################################################
### Set variables
export REFERENCE=~/space/s1/lin.yuan/puma/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
FASTA=$REFERENCE
#Regions to exclude:
#export REPEATMASK=~/project-klohmuel/Bricei/reference/GCA_028023285.1_mBalRic1.hap2_genomic.fna_rm4.1.4_TRF_merged.bed
################################################################################
#Activate conda environment

#GATK v4.2.2.0
#HTSJDK Version: 2.24.1
#Picard Version: 2.25.4
#bedtools v2.27.1
#bwa0.7.17-r1188
#samtools 1.19
#htslib 1.19
conda activate gatk_env
################################################################################
#Generate indexes
### Generate index files and intervals, sizes lists
samtools faidx ${FASTA}
cut -f1 ${FASTA}.fai > ${FASTA}.list
cut -f1-2 ${FASTA}.fai  > ${FASTA}.sizes

java -jar $PICARD CreateSequenceDictionary REFERENCE=${FASTA} OUTPUT=${FASTA%.fna*}.dict 
bwa index -a bwtsw ${FASTA} 

#Generate chunks
### Make 25 Mb intervals for splitting jobs
SIZES=${FASTA}.sizes
#generate windows
bedtools makewindows -g ${SIZES} -w 25000000 -s 25000000  > ${FASTA}.intervals_25Mb.bed
#Split by chromosome or scaffold
split -l1 --numeric-suffixes=1 -a 3 ${FASTA}.intervals_25Mb.bed ${FASTA}.intervals_25Mb_
for f in ${FASTA}.intervals_25Mb_* ; do mv ${f} ${f}.bed ; done

#Count intervals
N=$(ls ${FASTA}.intervals_25Mb_*.bed | wc -l)
N=$((N+1))

mkdir intervals
mv ${FASTA}.intervals_25Mb_* intervals

###DepthOfCoverage

export NUMTHREADS=8

# -mbq : minimum Phred quality score
# -mmq : minimum mapping quality
# -omitBaseOutput : do not output per base details (to omit unecessary output)
# -omitIntervals : exclude intervals specified by file (to omit unecessary output)
# -o : prefix outtput 
# -rf MappingQualityUnavailable : exclude if there is no mapping quality available
# There are several ways to estimate coverage, this tool matches the haplotype caller algorithm


for bam in *bam 
do 
NAME=${bam%%.bam}
gatk DepthOfCoverage \
-RF GoodCigarReadFilter \
-RF NonZeroReferenceLengthAlignmentReadFilter \
-RF PassesVendorQualityCheckReadFilter \
-RF MappingQualityAvailableReadFilter \
-RF MappingQualityReadFilter \
--minimum-mapping-quality 30 \
--min-base-quality 20 \
--omit-depth-output-at-each-base \
-R ${REFERENCE} \
-L ${REFERENCE}.list \
-I ${NAME}.yag.bam \
-O ${NAME}yag.bam
done

### Generate gVCF files (per chromosome)
export NUMTHREADS=4


#!/bin/bash
for interval in intervals.bed
do
REGION=$(ls $(dirname ${REFERENCE})/intervals/*_${IDX}.bed)
gatk HaplotypeCaller \
-ERC BP_RESOLUTION \
--minimum-mapping-quality 30 \
--min-base-quality-score 25 \
-R ${REFERENCE} \
-L ${REGION} \
-I ${NAME}.yag.bam \
-O ${NAME}_${IDX}.g.vcf.gz

################################################################################
### JOINT VCF FILE PROCESSING
### CombineGVCFs
#Combine individuals into a single VCF

REGION=$(ls $(dirname ${REFERENCE})/intervals/*_${IDX}.bed)
gatk CombineGVCFs \
-R ${REFERENCE} \
-L ${REGION} \
$(for i in *_${IDX}.g.vcf.gz ; do echo "-V ${i} "; done) \
-O bricei_joint_${IDX}.g.vcf.gz

### GenotypeGVCFs
gatk GenotypeGVCFs \
-all-sites \
-stand-call-conf 0 \
--only-output-calls-starting-in-intervals \
-R ${REFERENCE} \
-L ${REGION} \
-V bricei_joint_${IDX}.g.vcf.gz \
-O bricei_joint_${IDX}.vcf.gz

### LeftAlignAndTrimVariants
​
#Left Alignment: Variants in VCF files can sometimes have different representations, 
#especially when dealing with complex indels (insertions or deletions) or multiple 
#variant alleles. Left alignment standardizes the representation of variants by shifting
#them to the left-most position. This ensures that equivalent variants are represented consistently.

#Variant Trimming: Variant trimming removes any unnecessary reference bases on the left
#side of an indel variant. This can result in more concise and easily interpretable 
#variant representations in the VCF file
IDX=$(printf %03d ${SGE_TASK_ID})
REGION=$(ls $(dirname ${REFERENCE})/intervals/*_${IDX}.bed)
VCF=bricei_joint_${IDX}.vcf.gz
gatk LeftAlignAndTrimVariants \
-R ${REFERENCE} \
-L ${REGION} \
-V ${VCF} \
-O ${VCF%.vcf.gz}_LeftAlignTrim.vcf.gz


### SnpEff annotation v5.2
# Note: Following instructions in SnpEff manual, custom SnpEff annotation database created
# Jacqueline created the database
#Set up before using
echo "mBalRic1.hap2.genome : Balaenoptera_ricei" >> snpEff.config
echo "mBalRic1.hap2.reference : GCF_028023285.1_mBalRic1.hap2" >> snpEff.config
cp /u/project/klohmuel/DataRepository/Rice_whale/snpEff_mBalRic1.hap2.tar.gz ~/project-klohmuel/programs/snpEff/data
tar -xzvf ~/project-klohmuel/programs/snpEff/data/*.tar.gz

#Run snpEff
module load java/jdk-11.0.14
module load htslib 

IDX=$(printf %03d ${SGE_TASK_ID})

export SNPEFF=~/project-klohmuel/programs/snpEff/snpEff.jar
export DATABASE=mBalRic1.hap2
export VCF=bricei_joint_${IDX}_LeftAlignTrim.vcf.gz
export OUT=${VCF%_Left*}_snpEff.vcf.gz
export CSV=${OUT%.vcf.gz}_summary.csv

java -jar ${SNPEFF} -v -nodownload -csvStats ${CSV} ${DATABASE} ${VCF} \
| bgzip > ${OUT}

tabix -p vcf ${OUT}

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
