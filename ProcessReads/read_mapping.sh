ref=/space/s1/lin.yuan/puma/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
for sample in /space/s1/lin.yuan/puma/trimmed_reads/*.R1*
do name=${sample%.trim*}
var=$(echo $name |cut -d "/" -f7)
output=${var%_S*}
echo $output
bwa mem -t 15 $ref $name*R1* $name*R2*>$output.sam
samtools view -Sb -F 1804 $output.sam|samtools sort -@ 15 - >$output.bam
samtools index $output.bam
rm $output.sam
done



for bam in *bam ; do mv $bam ${bam%%.bam}.yag.bam;done
for bam in *bai ; do mv $bam ${bam%%.bam.bai}.yag.bam.bai;done
