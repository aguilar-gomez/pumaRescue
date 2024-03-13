#Add sample name to file
for file in *.bam
do
  name=$(basename $file .yag.bam);
  samtools addreplacerg -r "SM:${name}" -r "ID:xxxx" -o ${name}.y.bam $file -@ 10 ;
done

