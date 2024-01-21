for sample in *fastq.gz*
do
name=${sample%_*1*}
echo $name $sample ${name}_2.fastq.gz
java -jar ~/bin/trimmomatic-0.39.jar PE -threads 1 $sample ${name}_2.fastq.gz $name.trim.R1.fq.gz output_forward_unpaired.fq.gz $name.trim.R2.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:NexteraAdapters.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:75
rm *unpaired*
done


for sample in SRR13887964.trim.R1.fq
do
name=${sample%trim*1*}
echo $name $sample ${name}trim.R2.fq
perl /home/diana/bin/prinseq-lite.pl -fastq2 $sample ${name}trim.R2.fq -custom_params “AACCCT 5; G 15” -stats_all -out_good 
done
