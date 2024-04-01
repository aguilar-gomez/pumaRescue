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

#Tutorial:
#https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
#/u/home/d/daguilar/.conda/envs/repeatmask/bin/RepeatMasker - 4.1.5
#Mammalia runs well
RepeatMasker -engine ncbi -s -align -species "mammalia" -dir PyagRep $FASTA 

#Test runs well
RepeatMasker -engine ncbi -s -align -species "carnivora" -dir PyagRep $FASTA 
RepeatMasker -engine ncbi -s -align -species "felidae" -dir PyagRepFelidae $FASTA 



