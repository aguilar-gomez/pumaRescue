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

##############################################################################
### JOINT VCF FILE PROCESSING
#Combine individuals into a single VCF
#!/bin/sh

conda activate gatk_env
REFERENCE=/space/s1/lin.yuan/puma/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna

gatk CombineGVCFs \
  -R ${REFERENCE} \
  $(for i in *.vcf.gz ; do echo "-V ${i} "; done) \
  -O puma_allSamples.g.vcf.gz

### GenotypeGVCFs
gatk GenotypeGVCFs \
  -all-sites \
  -stand-call-conf 0 \
  --only-output-calls-starting-in-intervals \
  -R ${REFERENCE} \
  -V puma_allSamples.g.vcf.gz \
  -O puma_allSamples.vcf.gz
  
### LeftAlignAndTrimVariants
​
#Left Alignment: Variants in VCF files can sometimes have different representations, 
#especially when dealing with complex indels (insertions or deletions) or multiple 
#variant alleles. Left alignment standardizes the representation of variants by shifting
#them to the left-most position. This ensures that equivalent variants are represented consistently.

#Variant Trimming: Variant trimming removes any unnecessary reference bases on the left
#side of an indel variant. This can result in more concise and easily interpretable 
#variant representations in the VCF file
gatk LeftAlignAndTrimVariants \
  -R ${REFERENCE} \
  -V puma_allSamples.vcf.gz \
  -O puma_allSamples_LeftAlignTrim.vcf.gz

