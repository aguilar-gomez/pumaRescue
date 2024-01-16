for sample in *SRR7690239*
do
name=${sample%_*1*}
echo $name $sample ${name}_2.fastq.gz
java -jar ~/bin/trimmomatic-0.39.jar PE -threads 1 $sample ${name}_2.fastq.gz $n
ame.trim.R1.fq.gz output_forward_unpaired.fq.gz $name.trim.R2.fq.gz output_rever
se_unpaired.fq.gz ILLUMINACLIP:NexteraAdapters.fa:2:30:10:2:keepBothReads LEADIN
G:3 TRAILING:3 MINLEN:75
rm *unpaired*
done
