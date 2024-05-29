awk '$3 == "gene" {print $1, $4, $5, $7, $9}' OFS="\t" genomic.gff > geneCoordinates.tsv

export gff=/u/home/d/daguilar/project-kirk-bigdata/Pconcolor/genome_outgroup/annotation/genomic.gff
awk '$3 == "exon" {print $1, $4, $5, $7, $9}' OFS="\t" $gff > exonCoordinates.tsv

#In scaffolds that are at least 100kb
#Make windows 100kb, sliding 10kb
module load bedtools 
bedtools makewindows -g scaffolds100Kb.len -w 100000 -s 10000 > scaffolds_w100kb_s10kb
bedtools coverage -a scaffolds_w100kb_s10kb -b geneCoordinates.tsv  | awk '$6==100000' > Scaffold_window_GeneDensity.bed
bedtools intersect -b Scaffold_window_GeneDensity.bed -a TXprop.bed -wb|cut -f1-4,11|uniq > TXprop_geneD.bed
bedtools intersect -b Scaffold_window_GeneDensity.bed -a TXprop.bed -wb|cut -f1-4,8|uniq > TXprop_geneC.bed

#5mb windows
module load bedtools 
bedtools makewindows -g scaffolds100Kb.len -w 5000000 -s 1000000 > scaffolds_w5Mb_s1Mb
bedtools coverage -a scaffolds_w5Mb_s1Mb -b geneCoordinates.tsv  | awk '$6==5000000' > Scaffold_window_GeneDensity_w5Mb_s1Mb.bed
#find how much gene density is around each of the SNPs with called ancestry
bedtools intersect -b Scaffold_window_GeneDensity_w5Mb_s1Mb.bed -a TXprop.bed -wb|cut -f1-4,11|uniq > TXprop_geneD.bed
bedtools intersect -b Scaffold_window_GeneDensity_w5Mb_s1Mb.bed -a TXprop.bed -wb|cut -f1-4,8|uniq > TXprop_geneC.bed
#average ancestry within each window
bedtools intersect -a Scaffold_window_GeneDensity_w5Mb_s1Mb.bed -b TXprop.bed -wb -wa|cut -f1-4,7,11 > Scaffold_w5Mb_s1Mb_TXprop.bed


#5Mb windows, no sliding and exons
#5mb windows
module load bedtools 
bedtools makewindows -g scaffolds100Kb.len -w 5000000 > scaffolds_w5Mb
bedtools coverage -a scaffolds_w5Mb -b exonCoordinates.tsv  | awk '$6==5000000' > Scaffold_window_ExonDensity_w5Mb.bed
#average ancestry within each window
bedtools intersect -a Scaffold_window_ExonDensity_w5Mb.bed -b TXprop.bed -wb -wa|cut -f1-4,7,11 > ScaffoldExonD_w5Mb_TXprop.bed

bedtools coverage -a scaffolds_w5Mb -b geneCoordinates.tsv  | awk '$6==5000000' > Scaffold_window_GeneDensity_w5Mb.bed
bedtools intersect -a Scaffold_window_GeneDensity_w5Mb.bed -b TXprop.bed -wb -wa|cut -f1-4,7,11 > ScaffoldGeneD_w5Mb_TXprop.bed

#1Mb windows, no sliding and gene
module load bedtools 
bedtools makewindows -g scaffolds100Kb.len -w 1000000 > scaffolds_w1Mb
bedtools coverage -a scaffolds_w1Mb -b geneCoordinates.tsv  | awk '$6==1000000' > Scaffold_window_Density_w1Mb.bed
#average ancestry within each window
bedtools intersect -a Scaffold_window_Density_w1Mb.bed -b TXprop.bed -wb -wa|cut -f1-4,7,11 > Scaffold_w1Mb_TXprop.bed
