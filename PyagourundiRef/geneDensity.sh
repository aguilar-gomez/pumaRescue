awk '$3 == "gene" {print $1, $4, $5, $7, $9}' OFS="\t" genomic.gff > geneCoordinates.tsv

#In scaffolds that are at least 100kb
#Make windows 100kb, sliding 10kb
module load bedtools 
bedtools makewindows -g scaffolds100Kb.len -w 100000 -s 10000 > scaffolds_w100kb_s10kb
bedtools coverage -a scaffolds_w100kb_s10kb -b geneCoordinates.tsv  | awk '$6==100000' > Scaffold_window_GeneDensity.bed

bedtools intersect -b Scaffold_window_GeneDensity.bed -a TXprop.bed -wb|cut -f1-4,11|uniq > TXprop_geneD.bed
bedtools intersect -b Scaffold_window_GeneDensity.bed -a TXprop.bed -wb|cut -f1-4,8|uniq > TXprop_geneC.bed


module load bedtools 
bedtools makewindows -g scaffolds100Kb.len -w 5000000 -s 1000000 > scaffolds_w5Mb_s1Mb
bedtools coverage -a scaffolds_w5Mb_s1Mb -b geneCoordinates.tsv  | awk '$6==5000000' > Scaffold_window_GeneDensity_w5Mb_s1Mb.bed

bedtools intersect -b Scaffold_window_GeneDensity_w5Mb_s1Mb.bed -a TXprop.bed -wb|cut -f1-4,11|uniq > TXprop_geneD.bed
bedtools intersect -b Scaffold_window_GeneDensity_w5Mb_s1Mb.bed -a TXprop.bed -wb|cut -f1-4,8|uniq > TXprop_geneC.bed
