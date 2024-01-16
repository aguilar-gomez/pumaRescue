for sample in *SRR7148347*
do
echo $sample
fastqc -o fastqc_result -f fastq -t 5 $sample
done
