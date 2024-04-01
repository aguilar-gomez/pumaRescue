#!/bin/bash
# Repeat masking
# Define variables:
export FASTA=
export SPECIES="Pyaguarundi"
export NAME=${FASTA%.fa*}

source ~/anaconda3/etc/profile.d/conda.sh 

module load anaconda3
#conda create -n repeatmask
conda activate repeatmask
#conda install -c bioconda -c conda-forge repeatmasker

#Tutorial:
#https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

#/u/home/d/daguilar/.conda/envs/repeatmask/bin/RepeatMasker - 4.1.5
RepeatMasker -engine crossmatch -s -align -species "${SPECIES}" -dir $NAME $FASTA




