# copy site files
cp -s /space/s1/lin.yuan/puma/analysis_yagAsReference/down_sampling/all2yag_aftermaf_ohana* .


# prune site
PRUNEPARAM="10kb 1 .16"
IPNUT_SITE=all2yag_aftermaf_ohana
OUTPUT_NAME=all2yag_aftermaf_ahmm  
plink_1.90_6.26 --file $IPNUT_SITE --indep-pairwise $PRUNEPARAM --threads 10 --out $OUTPUT_NAME --allow-extra-chr



(base) lin.yuan@tilden:/space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm_May2024_4CYP$ wc -l all2yag_aftermaf_ahmm.prune.in
2842316 all2yag_aftermaf_ahmm.prune.in
(base) lin.yuan@tilden:/space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm_May2024_4CYP$ wc -l all2yag_aftermaf_ahmm.prune.out
17122048 all2yag_aftermaf_ahmm.prune.out



# prepare bamlist
(base) lin.yuan@tilden:/space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm_May2024_4CYP$ cp /space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm/puma_yagasref.bamlist .


# Change plink format to tab-separated sites
 while read line; do modified=${line::-3}; echo $modified ; done < all2yag_aftermaf_ahmm.prune.in|tr ":" "\t" > pruned_sites_forvcf


# Generate sample list for bcftools query
grep "AFP\|CYP\|TX" puma_yagasref.bamlist |awk '{print $0}' > sample_list4ahmm


# prepare input vcf
cp -s /space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm/filtered_all2yag.vcf.gz* .



# subset vcf
(base) lin.yuan@tilden:/space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm_May2024_4CYP$ nano subset_vcf.sh

# Run bcftools query to get partial samples (TX/CYP/AFP) for A_HMM and format it
#Nine fields before samples, include header
SAMPLE_LIST=sample_list4ahmm
REGION_FILE=pruned_sites_forvcf
INPUT_VCF=filtered_all2yag.vcf.gz
OUTPUT_NAME=all2yag_pruned4ahmm.vcf
OUTPUT_NAME2=all2yag_pruned4ahmm_unformat.vcf


bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%ID[\t%GT:%AD]\n' -R $REGION_FILE -S $SAMPLE_LIST $INPUT_VCF -o $OUTPUT_NAME &

bcftools view -o $OUTPUT_NAME2 -R $REGION_FILE -S $SAMPLE_LIST $INPUT_VCF --threads 15 &


# The header of “all2yag_pruned4ahmm.vcf” has [int] and we don’t want that. We will manually modify it (delete “[integer]”).

grep "#CHROM" all2yag_pruned4ahmm_unformat.vcf > newheader
nano newheader
# manually delete all the spaces and rename col9

# then concatenate two things
tail -n+2 all2yag_pruned4ahmm.vcf > noheader.vcf 
cat newheader noheader.vcf > all2yag_pruned_reheaded.vcf &



# create “sample2population” file
# 4CYP, 5TX
(base) lin.yuan@tilden:/space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm_May2024_4CYP$ more sample2population
AFP10.yag.bam   admixed
AFP11.yag.bam   admixed
AFP12.yag.bam   admixed
AFP13.yag.bam   admixed
AFP14.yag.bam   admixed
AFP15.yag.bam   admixed
AFP16.yag.bam   admixed
AFP17.yag.bam   admixed
AFP18.yag.bam   admixed
AFP19.yag.bam   admixed
AFP1.yag.bam    admixed
AFP20.yag.bam   admixed
AFP21.yag.bam   admixed
AFP22.yag.bam   admixed
AFP23.yag.bam   admixed
AFP24.yag.bam   admixed
AFP25.yag.bam   admixed
AFP26.yag.bam   admixed
AFP27.yag.bam   admixed
AFP28.yag.bam   admixed
AFP29.yag.bam   admixed
AFP2.yag.bam    admixed
AFP30.yag.bam   admixed
AFP31.yag.bam   admixed
AFP3.yag.bam    admixed
AFP4.yag.bam    admixed
AFP5.yag.bam    admixed
AFP6.yag.bam    admixed
AFP7.yag.bam    admixed
AFP8.yag.bam    admixed
AFP9.yag.bam    admixed
CYP45.yag.bam   0
CYP47.yag.bam   0
CYP51.yag.bam   0
CYP60.yag.bam   0
TX101.yag.bam   1
TX105.yag.bam   1
TX106.yag.bam   1
TX107.yag.bam   1
TX108.yag.bam   1


running vcf2ahmm
# use the python script (my modified version; I debugged a typo) to convert to needed format
nohup python3 vcf2ahmm.py -v all2yag_pruned_reheaded.vcf -s sample2population --min_total 4 -m 0 --min_diff 0.2 > ahmm_input &


# arguments
# "-g 1" indicates that admixed sample genotypes (rather than allele counts) should be used
# "-r [float]" sets the per site generation rate
# "-m [int]" minimum distance in basepairs between successful SNPs to be included in this analysis
# "-o [string]" file to print admixed sample ploidy for ahmm input. Default is ploidy.txt
# "--min_total [int]" minimum number of samples in each ancestral population to consider a site. Default 10.
# "--min_diff [float]" minimum allele frequency difference between any pair of ancestral populations to include a site. I.e., this selects AIMs. Default 0.1.
# there is site reduction here, from 2842316 to 284228

(base) lin.yuan@tilden:/space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm_May2024_4CYP$ wc -l ahmm_input
284228 ahmm_input
