#Add sample name to file
for file in *.bam
do
  name=$(basename $file .yag.bam);
  samtools addreplacerg -r "SM:${name}" -r "ID:xxxx" -o ${name}.y.bam $file -@ 10 ;
done

##############################################################################
#!/bin/sh
conda init bash
conda activate gatk_env

REFERENCE=/space/s1/lin.yuan/puma/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna

for file in *.y.bam
do
  name=$(basename $file .y.bam);
  samtools index ${file}

  gatk HaplotypeCaller -ERC BP_RESOLUTION \
  --minimum-mapping-quality 30 \
  --min-base-quality-score 25 \
  -R ${REFERENCE} \
  -I $name.y.bam \
  -O ${name}.vcf.gz 

done
##############################################################################
