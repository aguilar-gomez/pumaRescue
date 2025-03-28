#!/bin/bash
# Repeat masking
# Define variables:
export FASTA=GCF_014898765.1_PumYag_genomic.fna
export NAME=${FASTA%.fa*}

module load anaconda3
#conda create -n repeatmask
conda activate repeatmask
#conda install -c bioconda -c conda-forge repeatmasker

##Extract fasta that we are actually using
cat ${FASTA}_large_scaffolds.sizes ${FASTA}_small_scaffolds.bed|cut -f1 > scaffolds2keep
xargs samtools faidx $FASTA < scaffolds2keep > $NAME.reduced.fasta

#Genome is softmasked
# FROM NCBI:
#Are repetitive sequences in eukaryotic genomes masked?
#Repetitive sequences in eukaryotic genome assembly sequence files, 
#as identified by WindowMasker, have been masked to lower-case.
#The location and identity of repeats found by RepeatMasker are also
#provided in a separate file. These spans could be used to mask the
#genomic sequences if desired. Be aware, however, that many less 
#studied organisms do not have good repeat libraries available for 
#RepeatMasker to use.
python3 fasta_regex.py GCF_014898765.1_PumYag_genomic.fna.reduced.fasta "[atgcn]+" GCF_014898765.1_PumYag_genomic.SoftMask.bed
 
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o RM.log.txt
#$ -l highp,h_rt=72:00:00,h_data=24G
## and the number of cores as needed:
#$ -pe shared 4
#$ -M daguilar
#Run with specific species
export FASTA=GCF_014898765.1_PumYag_genomic.fna
export NAME=${FASTA%.fa*}
. /u/local/Modules/default/init/modules.sh

module load anaconda3
#conda create -n repeatmask
conda activate repeatmask
RepeatMasker -engine ncbi -s -align -species "jaguarundi" -dir PyagRepeatMasked $NAME.reduced.fasta -pa 4

#Repeat Masker output to bedfile
awk  -F ' ' 'BEGIN{ OFS = "\t" }{
    if (NR < 4) {}
    else {$6 = $6 - 1 ; print $5,$6,$7}
}'  GCF_014898765.1_PumYag_genomic.fna.reduced.fasta.out  > GCF_014898765.1_PumYag_genomic.fna.RMasker.bed


#Remove sex chromosome windows:
#Sex chromosome windows provided by Lin Yuan (male vs female coverage & annotation)
#First index
#mv windows2remove_both_Feb10_noColName.bed  sexChromosome_windows.bed
qrsh -l h_rt=4:00:00,h_data=40G -pe shared 2
module load gatk/4.2.0.0
module load bedtools

bedtools sort -i sexChromosome_windows.bed> sexChromosome_windows_sorted.bed
SOFTMASK=~/project-kirk-bigdata/Pconcolor/genome_outgroup/Masks/GCF_014898765.1_PumYag_genomic.SoftMask.bed
REPEATMASK=~/project-kirk-bigdata/Pconcolor/genome_outgroup/Masks/GCF_014898765.1_PumYag_genomic.fna.RMasker.bed 
SEXMASK=~/project-kirk-bigdata/Pconcolor/genome_outgroup/Masks/sexChromosome_windows_sorted.bed
cat $SEXMASK $SOFTMASK $REPEATMASK |bedtools sort | bedtools merge > RepeatSoftSexMasks.bed 
gatk IndexFeatureFile -I RepeatSoftSexMasks.bed 
#Total sequence to mask: 1,283,646,252

### SexMask VCF
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
SEXMASK=~/project-kirk-bigdata/Pconcolor/genome_outgroup/Masks/RepeatSoftSexMasks.bed 
VCF=puma_allsamples_${IDX}_snpEff.vcf.gz


gatk VariantFiltration \
-R ${REFERENCE} \
-L ${REGION} \
-mask ${SEXMASK} --mask-name "FAIL_Repeat_Sex" \
-verbosity ERROR \
-V $VCF \
-O ${VCF%.vcf.gz}_Mask_noSex.vcf.gz
