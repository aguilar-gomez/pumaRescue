genome=/space/s1/lin.yuan/puma/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
name=yag_allsites
# samtools faidx $genome
# no need to faidx cuz it is already indexed
cut -f1-2 $genome.fai >$name.length
bedtools makewindows -g $name.length -w 100000 -s 10000 >${name}_windows.bed
 
for name in *.bam
do
echo $name
samtools bedcov /space/s1/lin.yuan/puma/analysis/puma_coverage_calculation_map2yagGenome/yag_allsites_windows.bed $name > ${name}_window
done
