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
