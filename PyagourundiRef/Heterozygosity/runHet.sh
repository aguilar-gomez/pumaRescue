VCF=unfiltered_all2yag.vcf.gz

tabix -p vcf $VCF

while read scaffold len
do
echo $scaffold
bcftools view -r $scaffold -Oz -o pumaUnfil_${scaffold}.vcf.gz $VCF
tabix -p puma_${scaffold}.vcf.gz 
python3 ./SlidingWindowHet.py pumaUnfil_${scaffold}.vcf.gz scaffolds_10Mb 1000000 100000 $scaffold
done <scaffolds_10Mb
