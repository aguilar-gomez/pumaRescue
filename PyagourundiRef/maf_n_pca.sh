INPUT=all2yag_noMXY_windowfiltered.vcf
OUTPUT_NAME=all2yag_after_maf
plink_1.90_6.26 --vcf $INPUT --make-bed --out $OUTPUT_NAME --threads 16 --const-fid 0 --set-missing-var-ids @:#\$1,\$2 --snps-only --biallelic-only --allow-extra-chr --maf 0.01 --geno 0.0 --pca
