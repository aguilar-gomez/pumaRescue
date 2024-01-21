for sample in *bam
do
samtools depth  $sample| awk -v awkvar="$sample" 'BEGIN {OFS="\t"}; {sum+=$3; sumsq+=$3*$3} END { print awkvar, sum/NR, sqrt(sumsq/NR - (sum/NR)**2)}'
done
