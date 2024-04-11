#Previous
VCF=unfiltered_all2yag.vcf.gz

tabix -p vcf $VCF

while read scaffold len
do
echo $scaffold
bcftools view -r $scaffold -Oz -o pumaUnfil_${scaffold}.vcf.gz $VCF
tabix -p puma_${scaffold}.vcf.gz 
python3 ./SlidingWindowHet.py pumaUnfil_${scaffold}.vcf.gz scaffolds_10Mb 1000000 100000 $scaffold
done <scaffolds_10Mb

#################################################################################################################
#GATK
module load htslib
VCF=puma_simplePASS_all.vcf.gz
tabix -p vcf $VCF 
#Create file with scaffolds
intervals=~/project-kirk-bigdata/Pconcolor/genome_outgroup/intervals/*
cat $intervals|cut -f1|uniq > scaffolds_100Kb
cut -f1-2 GCF_014898765.1_PumYag_genomic.fna.reduced.fasta.fai > GCF_014898765.1_PumYag_genomic.fna.reduced.fasta.lengths
cp ~/project-kirk-bigdata/Pconcolor/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna.reduced.fasta.lengths scaffolds100Kb.len

#################################################################################################################
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o het.out
#$ -l highp,h_rt=72:00:00,h_data=24G
#$ -pe shared 2
#$ -M daguilar
#$ -t 1-176:1
. /u/local/Modules/default/init/modules.sh

module load bcftools
module load htslib 
VCF=puma_simplePASS_all.vcf.gz

awk -v var="${SGE_TASK_ID}" 'NR==var' scaffolds_100Kb > current_scaffold
read scaffold < current_scaffold
bcftools view -r $scaffold -Oz -o puma_${scaffold}.vcf.gz $VCF
tabix -p puma_${scaffold}.vcf.gz 
python3 ./SlidingWindowHet.py puma_${scaffold}.vcf.gz scaffolds100Kb.len 100000 10000 $scaffold

#################################################################################################################


