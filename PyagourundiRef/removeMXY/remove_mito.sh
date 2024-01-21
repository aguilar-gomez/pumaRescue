# remove by chromosome
# remove mito
# mito chr name got from here: https://www.ncbi.nlm.nih.gov/genome/?term=jaguarundi
# NC_028311.1
INPUT=/space/s1/lin.yuan/puma/analysis_yagAsReference/filtered_vcf/filtered_all2yag.vcf
OUTPUT_NAME=all2yag_noM
plink_1.90_6.26 --vcf $INPUT --make-bed --out $OUTPUT_NAME --threads 16 --set-missing-var-ids @:#\$1,\$2 --allow-extra-chr --not-chr NC_028311.1
