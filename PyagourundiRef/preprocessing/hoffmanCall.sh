### Bricei calling haplotypes and annotate(example commands)
# Workflow
# - Index genome and generate chunks
# - Per sample:
#   - DepthOfCoverage
#   - HaplotypeCaller
# - All samples jointly:
#   - GenotypeGVCFs 
#   - LeftAlignAndTrimVariants
#   - SnpEff
################################################################################
### Set variables
export REFERENCE=~/project-kirk-bigdata/Pconcolor/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
FASTA=$REFERENCE
#Regions to exclude:
#export REPEATMASK=CURRENTLY_DONT_HAVE_THIS
################################################################################
module load picard_tools/2.25.0 
module load samtools/1.15
module load bwa/0.7.17
module load bedtools/2.30.0
module load htslib
module load gatk/4.2.0.0
################################################################################
#Generate indexes
### Generate index files and intervals, sizes lists
samtools faidx ${FASTA}
cut -f1 ${FASTA}.fai > ${FASTA}.list
cut -f1-2 ${FASTA}.fai  > ${FASTA}.sizes

picard CreateSequenceDictionary REFERENCE=${FASTA} OUTPUT=${FASTA%.fna*}.dict 
bwa index -a bwtsw ${FASTA} 

#Generate chunks
### Make 5 Mb intervals for splitting jobs
SIZES=${FASTA}.sizes
#generate windows
# Generate windows for scaffolds larger than 5 Mb
awk '$2 > 5000000 {print $0}' ${SIZES} > ${FASTA}_large_scaffolds.sizes
#generate windows for large scaffolds
bedtools makewindows -g ${FASTA}_large_scaffolds.sizes -w 25000000 -s 25000000  > ${FASTA}.intervals_25Mb.bed
# Generate a single interval for scaffolds <=5 Mb and they have to be >100Kb
awk '$2 <= 5000000 && $2 > 100000 {printf "%s\t%d\t%d\n", $1, 0, $2}' "${SIZES}" >${FASTA}_small_scaffolds.bed

cat ${FASTA}.intervals_25Mb.bed ${FASTA}_small_scaffolds.bed > ${FASTA}_intervals.bed
#TOTAL SEQUENCE IN INTERVALS
#awk '{sum += $3 - $2} END {print sum}' GCF_014898765.1_PumYag_genomic.fna_intervals.bed
#2,415,546,177

# Split large scaffolds by chromosome or scaffold
split -l 1 --numeric-suffixes=1 -a 3 ${FASTA}_intervals.bed ${FASTA}_intervals_
for f in ${FASTA}_intervals_* ; do mv ${f} ${f}.bed ; done


###DepthOfCoverage
# -mbq : minimum Phred quality score
# -mmq : minimum mapping quality
# -omitBaseOutput : do not output per base details (to omit unecessary output)
# -omitIntervals : exclude intervals specified by file (to omit unecessary output)
# -o : prefix outtput 
# -rf MappingQualityUnavailable : exclude if there is no mapping quality available
# There are several ways to estimate coverage, this tool matches the haplotype caller algorithm



module load gatk/4.2.0.0
module load samtools/1.15

REFERENCE=~/project-kirk-bigdata/Pconcolor/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
for bam in *.y.bam; 
do
samtools index $bam
NAME=${bam%%.y.*}

gatk DepthOfCoverage \
  -RF GoodCigarReadFilter \
  -RF NonZeroReferenceLengthAlignmentReadFilter \
  -RF PassesVendorQualityCheckReadFilter \
  -RF MappingQualityAvailableReadFilter \
  -RF MappingQualityReadFilter \
  --minimum-mapping-quality 30 \
  --min-base-quality 25 \
  --omit-depth-output-at-each-base \
  -R ${REFERENCE} \
  -L ${REFERENCE}.list \
  -I $bam \
  -O ${NAME}.DepOfCoverage
done

################################################################################
### Generate gVCF files 
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o HapCaller1.log.$JOB_ID.$TASK_ID
#$ -l highp,h_rt=72:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 4
#$ -M daguilar
#$ -t 1-237:1

. /u/local/Modules/default/init/modules.sh

module load gatk/4.2.0.0

NAME=$1
BAM=$NAME.y.bam

REFERENCE=~/project-kirk-bigdata/Pconcolor/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
IDX=$(printf %03d ${SGE_TASK_ID})
REGION=$(ls $(dirname ${REFERENCE})/intervals/*_${IDX}.bed)

gatk HaplotypeCaller \
  -ERC BP_RESOLUTION \
  --minimum-mapping-quality 30 \
  --min-base-quality-score 25 \
  -R ${REFERENCE} \
  -L ${REGION} \
  -I $BAM \
  -O ${NAME}_${IDX}.g.vcf.gz 

################################################################################
### JOINT VCF FILE PROCESSING
### CombineGVCFs
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o Combine.log.$JOB_ID.$TASK_ID
#$ -l highp,h_rt=72:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 4
#$ -M daguilar
#$ -t 1-237:1

. /u/local/Modules/default/init/modules.sh

module load gatk/4.2.0.0
#Combine individuals into a single VCF
IDX=$(printf %03d ${SGE_TASK_ID})
REFERENCE=~/project-kirk-bigdata/Pconcolor/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
REGION=$(ls $(dirname ${REFERENCE})/intervals/*_${IDX}.bed)

gatk CombineGVCFs \
-R ${REFERENCE} \
-L ${REGION} \
$(for i in *_${IDX}.g.vcf.gz ; do echo "-V ${i} "; done) \
-O puma_${IDX}.g.vcf.gz

### GenotypeGVCFs
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o genotype.log.$JOB_ID.$TASK_ID
#$ -l highp,h_rt=72:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 4
#$ -M daguilar
#$ -t 1-237:1
. /u/local/Modules/default/init/modules.sh

module load gatk/4.2.0.0
IDX=$(printf %03d ${SGE_TASK_ID})
REFERENCE=~/project-kirk-bigdata/Pconcolor/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
REGION=$(ls $(dirname ${REFERENCE})/intervals/*_${IDX}.bed)

gatk GenotypeGVCFs \
  -all-sites \
  -stand-call-conf 0 \
  --only-output-calls-starting-in-intervals \
  -R ${REFERENCE} \
  -L ${REGION} \
  -V puma_${IDX}.g.vcf.gz \
  -O puma_allsamples_${IDX}.vcf.gz

### LeftAlignAndTrimVariants
​
#Left Alignment: Variants in VCF files can sometimes have different representations, 
#especially when dealing with complex indels (insertions or deletions) or multiple 
#variant alleles. Left alignment standardizes the representation of variants by shifting
#them to the left-most position. This ensures that equivalent variants are represented consistently.

#Variant Trimming: Variant trimming removes any unnecessary reference bases on the left
#side of an indel variant. This can result in more concise and easily interpretable 
#variant representations in the VCF file
### GenotypeGVCFs
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o genotype.log.$JOB_ID.$TASK_ID
#$ -l highp,h_rt=72:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 2
#$ -M daguilar
#$ -t 1-237:1
. /u/local/Modules/default/init/modules.sh

module load gatk/4.2.0.0

IDX=$(printf %03d ${SGE_TASK_ID})
REFERENCE=~/project-kirk-bigdata/Pconcolor/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
REGION=$(ls $(dirname ${REFERENCE})/intervals/*_${IDX}.bed)
VCF=puma_allsamples_${IDX}.vcf.gz

gatk LeftAlignAndTrimVariants \
  -R ${REFERENCE} \
  -L ${REGION} \
  -V ${VCF} \
  -O ${VCF%.vcf.gz}_LeftAlignTrim.vcf.gz


### SnpEff annotation 
# Note: Following instructions in SnpEff manual, custom SnpEff annotation database created
#Set up before using
#SNPeff database
#on tilden
#DATABASE=/space/s1/lin.yuan/puma/analysis_yagAsReference/snpeff/data/PumYag
tar -czvf PumYag.tar.gz PumYag

#hoffman
cd ~/project-klohmuel/programs/snpEff
echo "PumYag.genome : Puma_yaguarundi" >> snpEff.config
echo "PumYag.reference : GCF_014898765.1_PumYag_genomic" >> snpEff.config
cp ~/project-kirk-bigdata/Pconcolor/PumYag.tar.gz ~/project-klohmuel/programs/snpEff/data
tar -xzvf ~/project-klohmuel/programs/snpEff/data/PumYag.tar.gz

#Run snpEff

#!/bin/bash
#$ -cwd
#$ -j y
#$ -o snpEFF.log.$JOB_ID.$TASK_ID
#$ -l highp,h_rt=24:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 2
#$ -M daguilar
#$ -t 1-237:1
. /u/local/Modules/default/init/modules.sh

module load java/jdk-11.0.14
module load htslib 

IDX=$(printf %03d ${SGE_TASK_ID})

export SNPEFF=~/project-klohmuel/programs/snpEff/snpEff.jar
export DATABASE=PumYag
export VCF=puma_allsamples_${IDX}_LeftAlignTrim.vcf.gz
export OUT=${VCF%_Left*}_snpEff.vcf.gz
export CSV=${OUT%.vcf.gz}_summary.csv

java -jar ${SNPEFF} -v -nodownload -csvStats ${CSV} ${DATABASE} ${VCF} \
| bgzip > ${OUT}

tabix -p vcf ${OUT}
