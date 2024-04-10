### Runs of homozygosity (ROH)
#!/bin/sh
#$ -l highp,h_rt=72:00:00,h_data=24G
#$ -pe shared 4
#$ -N ROHs
#$ -cwd
#$ -m bea
#$ -o ./rohv2.out
#$ -e ./rohv2.err
#$ -M daguilar
VCF=puma_simplePASS_variants_all.vcf.gz 
bcftools query -l $VCF >samples.list

# Identify ROH
# -G, --GTs-only <float>             use GTs and ignore PLs, instead using <float> for PL of the two least likely genotypes.
bcftools roh -G 30 -Orz -o ${VCF}_bcftoolsROH.txt.gz ${VCF}

#Repeat excluding related individuals and using only autosomes
#KEEP=samples_keep.list # generate based on relatedness
# Identify ROH
#bcftools roh -e ${KEEP} -G 30 -Orz -o ${VCF}_bcftoolsROH.txt.gz ${VCF}

# Reformat output
zcat ${VCF}_bcftoolsROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> ${VCF}_bcftoolsROH.txt

#NUMBER OF VARIANT AND INVARIANT SITES THAT PASS FILTERS 1059718553 THIS NUMBER NEEDS TO BE MODIFIED
# Calculate Froh using max length calculated from pseudo-genome
DATA=${VCF}_bcftoolsROH.txt.gz
while read -r SAMPLE ; do 
zcat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e6){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/1059718553, num, sum/num}'
done < samples.list

#!/bin/sh
#$ -l highp,h_rt=72:00:00,h_data=24G
#$ -pe shared 4
#$ -N ROHs
#$ -cwd
#$ -m bea
#$ -o ./rohv2.out
#$ -e ./rohv2.err
#$ -M daguilar

. /u/local/Modules/default/init/modules.sh

module load bcftools
module load htslib


VCF=bricei_simplePASS_variants_autosomes.vcf.gz
#zcat ${VCF} | head -n 1000 | grep "^#" | tail -n 1 | cut -f10- | tr '\t' '\n' > samples.list
#KEEP=samples_keep.list
#grep -v "Bricei011\| Bricei016\| Bricei003\| Bricei024\| Bricei001\| Bricei007\|Bricei021\| Bricei005\| Bricei042\| Bricei015\| Bricei023" samples.list > ${KEEP}
#grep -v "Bede011\|Bede016\|Bede003\|Bede024\|Bede001\|Bede007\|Bede021\|Bede005\|Bede042\|Bede015\|Bede023" samples.list > ${KEEP}
# Identify ROH
bcftools roh -e ${KEEP} -G 30 -Orz -o ${VCF}_bcftoolsROH.txt.gz ${VCF}

# Identify ROH
# -G, --GTs-only <float>             use GTs and ignore PLs, instead using <float> for PL of the two least likely genotypes.
bcftools roh -G 30 -Orz -o ${VCF}_bcftoolsROH.txt.gz ${VCF}
echo STEP1

# Reformat output
zcat ${VCF}_bcftoolsROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> ${VCF}_bcftoolsROH.txt

echo STEP2
# Run with a pseudo-genome with entirely 0/0 genotypes
#zcat ${IN} | head -n 1000 | grep "^#" > head.tmp
# Manually edit head.tmp to add a sample name to the last line ("pseudohom")
IN=bricei_simplePASS_variants_autosomes.vcf.gz
OUT=bricei_simplePASS_variants_autosomes_pseudohom.vcf.gz
zcat ${IN} | grep -v "^#" | sed -e 's/$/\t0\/0/g' | cat head.tmp - | bgzip > ${OUT}
bcftools roh -e ${KEEP} -G 30 -Orz -o bricei_simplePASS_variants_autosomes_pseudohom.vcf.gz ${OUT}
zcat ${OUT}_bcftoolsROH.txt.gz \
| awk -v s=pseudohom 'BEGIN{sum=0}{if ($2==s){sum+=$6}}END{printf "%s\t%s\n", s, sum}'


# Calculate Froh using max length calculated from pseudo-genome
DATA=${OUT}_bcftoolsROH.txt.gz
while read -r SAMPLE ; do 
zcat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e6){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/
2560651057, num, sum/num}'
done < samples.list


### vcftools see if you can run without splitting by chromosomes
VCF=bricei_simplePASS_variants_chr.vcf.gz
for CHR in chr{1..21} ; do
vcftools --gzvcf ${VCF} --LROH --chr ${CHR} --out ${VCF}_vcftools_${CHR}
done

OUT=${VCF}_vcftools.LROH
cp ${VCF}_vcftools_chr1.LROH  ${OUT}
cat ${VCF}_vcftools_chr{2..21}.LROH | grep -v "^CHROM" >> ${OUT}

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






