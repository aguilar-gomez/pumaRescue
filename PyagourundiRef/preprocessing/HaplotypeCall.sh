#Add sample name to file
for file in *.bam; do
base_name=$(basename $file .yag.bam);
samtools addreplacerg -r "SM:${base_name}" -r "ID:xxxx" -o ${base_name}.y.bam $file -@ 10 ;
done


conda activate gatk_env

REFERENCE=/space/s1/lin.yuan/puma/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna

for file in *.y.bam; do
name=$(basename $file .y.bam);
samtools index ${name}

gatk HaplotypeCaller -ERC BP_RESOLUTION \
--minimum-mapping-quality 30 \
--min-base-quality-score 25 \
-R ${REFERENCE} \
-I $NAME.y.bam \
-O ${NAME}.vcf.gz 

done
