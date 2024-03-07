### Runs of homozygosity (ROH)
#Code from the vaquita project jarobin

bcftools query -l bricei_simplePASS_variants_all.vcf.gz >samples.list

module load bcftools
module load htslib
# Generate list of samples, exclude related individuals for calculating allele frequencies
# in bcftools roh command later
VCF=bricei_simplePASS_variants_all.vcf.gz
zcat ${VCF} | head -n 1000 | grep "^#" | tail -n 1 | cut -f10- | tr '\t' '\n' > samples.list
# Identify ROH
# -G, --GTs-only <float>             use GTs and ignore PLs, instead using <float> for PL of the two least likely genotypes.
bcftools roh -G 30 -Orz -o ${VCF}_bcftoolsROH.txt.gz ${VCF}

#Repeat excluding related individuals and using only autosomes
VCF=bricei_simplePASS_variants_autosomes.vcf.gz
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

#Bede001	0.555961	243	5.85853e+06
#Bede002	0.371319	270	3.52155e+06
#Bede003	0.431451	237	4.66158e+06
#Bede004	0.345505	235	3.76476e+06
#Bede005	0.375494	270	3.56114e+06
#Bede007	0.410193	235	4.46962e+06
#Bede008	0.384636	264	3.73075e+06
#Bede011	0.389005	254	3.92167e+06
#Bede015	0.411716	263	4.00859e+06
#Bede016	0.435883	253	4.41163e+06
#Bede021	0.386178	262	3.7743e+06
#Bede022	0.393442	262	3.8453e+06
#Bede023	0.405665	249	4.17175e+06
#Bede024	0.443978	254	4.47587e+06
#Bede034	0.371651	267	3.5643e+06
#Bede035	0.406754	269	3.87195e+06
#Bede037	0.386977	269	3.68369e+06
#Bede039	0.418523	280	3.82747e+06
#Bede040	0.342677	243	3.61102e+06
#Bede042	0.362281	255	3.63794e+06
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
#sample    totalROH    #segments
#Bede001 	 1161893190 	 864
#Bede002 	 822734054 	 966
#Bede003 	 957311132 	 836
#Bede004 	 780093483 	 891
#Bede005 	 807271928 	 974
#Bede007 	 890383648 	 757
#Bede008 	 896782323 	 855
#Bede011 	 825187204 	 800
#Bede015 	 916418677 	 980
#Bede016 	 900987350 	 833
#Bede021 	 859040983 	 771
#Bede022 	 830732384 	 823
#Bede023 	 900333867 	 755
#Bede024 	 968846251 	 755
#Bede034 	 837508895 	 891
#Bede035 	 845109632 	 790
#Bede037 	 856554476 	 828
#Bede039 	 919852541 	 917
#Bede040 	 776789174 	 749
#Bede042 	 782756525 	 863

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




