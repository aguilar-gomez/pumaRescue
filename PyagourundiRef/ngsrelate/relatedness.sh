#!/bin/sh
#$ -l highp,h_rt=72:00:00,h_data=24G
#$ -pe shared 4
#$ -N Relatedness
#$ -cwd
#$ -m bea
#$ -o ./relate.out
#$ -e ./relate.err
#$ -M daguilar

NGSRELATE=~/bin/ngsRelate/ngsRelate
$NGSRELATE -h puma_simplePASS_variants_all.vcf.gz -O relatedness_puma -T GT
