VCF=puma_simplePASS_variants_all.vcf.gz
cp -s ~/project-klohmuel/Pconcolor/vcfs/PASS8/$VCF .
module load plink
#PLINK v1.90b6.24
plink --vcf $VCF--make-bed --out puma_variants \
                --threads 15 --const-fid 0 --set-missing-var-ids @:#\$1,\$2 --snps-only \
                --biallelic-only --allow-extra-chr --vcf-half-call r

