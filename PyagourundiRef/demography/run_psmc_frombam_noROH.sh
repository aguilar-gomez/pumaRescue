#!/bin/bash

bamfile=$1
avdep=$2
bedfile=$3
mindep=$((avdep/3))
maxdep=$((2*avdep))
echo $mindep
echo $maxdep

# We filtered
# the alignment file for each puma to include sites which had between one third and twice
# the average coverage for that puma. We used a generation time of 5 years and a per
# generation mutation rate of 0.5e-8 per bp per generation

# For the fq2psmcfa step we set s =100 to account for the low heterozygosity of pumas
# We used the same parameters for the psmc step as used for human data
#${PSMC} -N25 -t15 -r5 -p "4+25*2+4+6"

out=$(basename $bamfile .bam)
{ time samtools mpileup -C50 -u -v -f data/puma.v4.1x.final.mt.fasta --positions $bedfile $bamfile 2> fq/${out}_stderr.txt  | bcftools call -c - | vcfutils.pl vcf2fq -d $mindep -D $maxdep | gzip > fq/${out}noROH.fq.gz ; } 2> fq/bam2fq_${out}noROH.time 
echo "Waiting for bam2fq to finish"
wait
echo "Running fq2psmcfa"
{ time fq2psmcfa -q 30 -s 100 fq/${out}noROH.fq.gz > psmcIN/${out}noROH.psmcfa ; } 2> psmcIN/fq2psmc_${out}noROH.time
{ time psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o psmcOUT/${out}noROH.psmc psmcIN/${out}noROH.psmcfa ; } 2> psmcOUT/psmc_${out}noROH.time
wait


# our code
# with ROH, no bootstrap
# command example: nohup sh psmc.sh TX101.yag.bam 11 NONE > out.pmsc_TX101 &
conda create -n samtools1.9
conda install -c bioconda samtools=1.9

out=$(basename $bamfile .bam)
{ time samtools mpileup -C50 -u -v -f GCF_014898765.1_PumYag_genomic.fna $bamfile 2> fq/${out}_stderr.txt  | bcftools call -c - | vcfutils.pl vcf2fq -d $mindep -D $maxdep | gzip > fq/${out}withROH.fq.gz ; } 2> fq/bam2fq_${out}withROH.time 
echo "Waiting for bam2fq to finish"
wait
echo "Running fq2psmcfa"
{ time /home/lin.yuan/programs/psmc-master/utils/fq2psmcfa -q 30 -s 100 fq/${out}withROH.fq.gz > psmcIN/${out}withROH.psmcfa ; } 2> psmcIN/fq2psmc_${out}withROH.time
{ time psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o psmcOUT/${out}withROH.psmc psmcIN/${out}withROH.psmcfa ; } 2> psmcOUT/psmc_${out}withROH.time
wait

