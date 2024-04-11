qrsh -l h_rt=4:00:00,h_data=40G -pe shared 2

VCF=puma_simplePASS_variants_all.vcf.gz
cp -s ~/project-klohmuel/Pconcolor/vcfs/PASS8/$VCF .

module load plink
#PLINK v1.90b6.24

#All variants no Maf
plink --vcf $VCF--make-bed --out puma_variants \
                --threads 15 --const-fid 0 --set-missing-var-ids @:#\$1,\$2 --snps-only \
                --biallelic-only --allow-extra-chr --vcf-half-call r

#MAF filter and PCA
OUTPUT_NAME=puma_variantsmaf1
plink --vcf $VCF --make-bed --out $OUTPUT_NAME --threads 16 --const-fid 0 \
                   --set-missing-var-ids @:#\$1,\$2 --snps-only --biallelic-only \
                   --allow-extra-chr --maf 0.01 --geno 0.05 --pca



OUTPUT_NAME=puma_variantsmaf5
plink --vcf $VCF --make-bed --out $OUTPUT_NAME --threads 16 --const-fid 0 \
                   --set-missing-var-ids @:#\$1,\$2 --snps-only --biallelic-only \
                   --allow-extra-chr --maf 0.05 --geno 0.0 --pca
#7,268,511 SNPs
#4,576,792 variants and 50 people pass filters and QC.

module load plink
#FST
awk '{print $1 "\t" $2 "\t" substr($2, 1, 2)}' puma_variantsmaf1.fam > puma.clust
awk 'BEGIN {OFS="\t"} {if ($3 == "AF") $3="FL"; if ($3 == "CY") $3="FL"; print}' puma.clust > puma.clustFL
awk 'BEGIN {OFS="\t"} {if ($3 == "SM") $3="CA"; if ($3 == "SC") $3="CA"; print}' puma.clustFL > puma.clustFLCA
#VCF=puma_simplePASS_variants_all.vcf.gz
plink2 --bfile puma_variantsmaf1 --fst "CATPHENO" 'method=wc'  \
        --within puma.clust --out puma --threads 10 --allow-extra-chr

plink2 --bfile puma_variantsmaf1 --fst "CATPHENO" 'method=wc'  \
        --within puma.clustFL --out pumaFL --threads 10 --allow-extra-chr

plink2 --bfile puma_variantsmaf1 --fst "CATPHENO" 'method=wc'  \
        --within puma.clustFLCA --out pumaFLCA --threads 10 --allow-extra-chr
