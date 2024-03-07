### Runs of homozygosity (ROH)
#modified code from jarobin

VCF=
bcftools query -l $VCF >samples.list

# Identify ROH
# -G, --GTs-only <float>             use GTs and ignore PLs, instead using <float> for PL of the two least likely genotypes.
bcftools roh -G 30 -Orz -o ${VCF}_bcftoolsROH.txt.gz ${VCF}

#Repeat excluding related individuals and using only autosomes
KEEP=samples_keep.list
grep -v "Bede011\|Bede016\|Bede003\|Bede024\|Bede001\|Bede007\|Bede021\|Bede005\|Bede042\|Bede015\|Bede023" samples.list > ${KEEP}
# Identify ROH
bcftools roh -e ${KEEP} -G 30 -Orz -o ${VCF}_bcftoolsROH.txt.gz ${VCF}

# Reformat output
zcat ${VCF}_bcftoolsROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> ${VCF}_bcftoolsROH.txt

# Run with a pseudo-genome with entirely 0/0 genotypes
IN=bricei_simplePASS_variants_all.vcf.gz
OUT=bricei_simplePASS_variants_pseudohom.vcf.gz

zcat ${IN} | head -n 1000 | grep "^#" > head.tmp
# Manually edit head.tmp to add a sample name to the last line ("pseudohom")
zcat ${IN} | grep -v "^#" | sed -e 's/$/\t0\/0/g' | cat head.tmp - | bgzip > ${OUT}
bcftools roh -G 30 -Orz -o ${OUT}_bcftoolsROH.txt.gz ${OUT}
zcat ${OUT}_bcftoolsROH.txt.gz | awk -v s=pseudohom 'BEGIN{sum=0}{if ($3==s){sum+=$6}}END{printf "%s\t%s\n", s, sum}'
zcat ${OUT}_bcftoolsROH.txt.gz | grep "pseudohom"|awk '{ sum += $6 } END { print sum }'

# Calculate Froh using max length calculated from pseudo-genome
DATA=${OUT}_bcftoolsROH.txt.gz
while read -r SAMPLE ; do 
zcat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e6){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/2560651057, num, sum/num}'
done < samples.list

#Do before
cat ${IN} | head -n 1000 | grep "^#" > head.tmp
# Manually edit head.tmp to add a sample name to the last line ("pseudohom")
#Autosomes and remove related
IN=bricei_simplePASS_variants_autosomes.vcf.gz
OUT=bricei_simplePASS_variants_autosomes_pseudohom.vcf.gz
zcat ${IN} | grep -v "^#" | sed -e 's/$/\t0\/0/g' | cat head.tmp - | bgzip > ${OUT}
bcftools roh -e ${KEEP} -G 30 -Orz -o ${OUT}_bcftoolsROH.txt.gz ${OUT}
zcat ${OUT}_bcftoolsROH.txt.gz \
| awk -v s=pseudohom 'BEGIN{sum=0}{if ($2==s){sum+=$6}}END{printf "%s\t%s\n", s, sum}'
#pseudohom	2,395,187,392
#pseudohom	2395187392

# Calculate Froh using max length calculated from pseudo-genome
DATA=${OUT}_bcftoolsROH.txt.gz
while read -r SAMPLE ; do 
zcat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e6){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/2395187392, num, sum/num}'
done < samples.list

### vcftools
VCF=bricei_simplePASS_variants_chr.vcf.gz
for CHR in chr{1..21} ; do
vcftools --gzvcf ${VCF} --LROH --chr ${CHR} --out ${VCF}_vcftools_${CHR}
done

OUT=${VCF}_vcftools.LROH
cp ${VCF}_vcftools_chr1.LROH  ${OUT}
cat ${VCF}_vcftools_chr{2..21}.LROH | grep -v "^CHROM" >> ${OUT}

grep "Bede" bricei_simplePASS_variants_chr.vcf.gz_vcftools.LROH |cut -f8|sort|uniq > samples.list
while read -r line
do
    cat bricei_simplePASS_variants_chr.vcf.gz_vcftools.LROH | grep "$line" | awk -v line="$line" '{ diff = $3 - $2; sum += diff } END { print line, "\t",sum, "\t", NR }'
done < samples.list

### plink
module load plink
# Convert simplified, concatenated VCF file to plink format
FILE=bricei_order
plink --vcf bricei_simplePASS_variants_all.vcf.gz --make-bed --out $FILE --allow-extra-chr --keep-allele-order --double-id 

# Identify ROH and reformat output
plink --bfile ${FILE} --out ${FILE}_ROH --homozyg --allow-extra-chr 
#--homozyg: Scan complete, found 5214 ROH.
sed -i -E 's/^\ +//g' ${FILE}_ROH.hom
sed -i -E 's/\ +/\t/g' ${FILE}_ROH.hom 

#KB       Length of region (kb)
#NSNP     Number of SNPs in run
#DENSITY  Average SNP density (1 SNP per kb)
#PHOM     Proportion of sites homozygous
#PHET     Proportion of sites heterozygous

#Autosomes
FILE=bricei_autosomes
plink --vcf bricei_simplePASS_variants_autosomes.vcf.gz --make-bed --out $FILE --chr-set 21 --keep-allele-order --double-id 
# Identify ROH and reformat output
plink --bfile ${FILE} --out ${FILE}_ROH --homozyg --allow-extra-chr 
#--homozyg: Scan complete, found 5214 ROH.
sed -i -E 's/^\ +//g' ${FILE}_ROH.hom
sed -i -E 's/\ +/\t/g' ${FILE}_ROH.hom 




