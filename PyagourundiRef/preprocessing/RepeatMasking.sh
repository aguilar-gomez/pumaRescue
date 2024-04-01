#!/bin/bash
# Repeat masking
# Define variables:
export FASTA=GCF_014898765.1_PumYag_genomic.fna
export SPECIES="Pyagouaroundi"
export NAME=${FASTA%.fa*}

source ~/anaconda3/etc/profile.d/conda.sh 

module load anaconda3
#conda create -n repeatmask
conda activate repeatmask
#conda install -c bioconda -c conda-forge repeatmasker

##Extract fasta that we are actually using
cat ${FASTA}_large_scaffolds.sizes ${FASTA}_small_scaffolds.bed|cut -f1 > scaffolds2keep
xargs samtools faidx $FASTA < scaffolds2keep > $NAME.reduced.fasta

#Genome is softmasked
python3 fasta_regex.py GCF_014898765.1_PumYag_genomic.fna.reduced.fasta "[atgcn]+" GCF_014898765.1_PumYag_genomic.SoftMask.bed


#Tutorial:
#https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
#/u/home/d/daguilar/.conda/envs/repeatmask/bin/RepeatMasker - 4.1.5
#Mammalia runs well
RepeatMasker -engine ncbi -s -align -species "mammalia" -dir PyagRep $FASTA 

#Run with specific species
RepeatMasker -engine ncbi -s -align -species "jaguarundi" -dir PyagRepeatMasked $NAME.reduced.fasta 

### RepeatMask VCF
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o RM.log.$JOB_ID.$TASK_ID
#$ -l highp,h_rt=72:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 1
#$ -M daguilar
#$ -t 1-237:1

. /u/local/Modules/default/init/modules.sh

module load gatk/4.2.0.0
module load htslib

REFERENCE=~/project-kirk-bigdata/Pconcolor/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
IDX=$(printf %03d ${SGE_TASK_ID})
REGION=$(ls $(dirname ${REFERENCE})/intervals/*_${IDX}.bed)

#Regions to exclude:
export REPEATMASK=~/project-kirk-bigdata/Pconcolor/genome_outgroup/Masks/GCF_014898765.1_PumYag_genomic.SoftMask.bed
### VariantFiltration to mask repeats
# Note: GATK requires indexed mask file
#gatk IndexFeatureFile -I ${REPEATMASK}

VCF=puma_allsamples_${IDX}_snpEff_filter.vcf.gz
tabix -p vcf $VCF

gatk VariantFiltration \
-R ${REFERENCE} \
-L ${REGION} \
-mask ${REPEATMASK} --mask-name "FAIL_Repeat" \
-verbosity ERROR \
-V $VCF \
-O ${VCF%.vcf.gz}_LeftAlignTrim_Mask.vcf.gz


