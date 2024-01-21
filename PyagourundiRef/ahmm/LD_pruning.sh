# copy site files
cp -s /space/s1/lin.yuan/puma/analysis_yagAsReference/down_sampling/all2yag_aftermaf_ohana* .




# prune site
PRUNEPARAM="10kb 1 .16"
IPNUT_SITE=all2yag_aftermaf_ohana
OUTPUT_NAME=all2yag_aftermaf_ahmm  
plink_1.90_6.26 --file $IPNUT_SITE --indep-pairwise $PRUNEPARAM --threads 20 --out $OUTPUT_NAME --allow-extra-chr


# --indep-pairwise <window size>['kb'] <step size (variant ct)> <r^2 threshold>
# keep same as our old runs


lin.yuan@tilden:/space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm$ wc all2yag_aftermaf_ahmm.prune.in
 2842316  2842316 75359189 all2yag_aftermaf_ahmm.prune.in
lin.yuan@tilden:/space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm$ wc all2yag_aftermaf_ahmm.prune.out
 17122048  17122048 455335250 all2yag_aftermaf_ahmm.prune.out


# prepare bamlist
cd /space/s1/lin.yuan/puma/bam_output_allsampletoOutgroup/allbams/
ls *.bam > puma_yagasref.bamlist
mv puma_yagasref.bamlist /space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm/puma_yagasref.bamlist
# Change plink format to tab-separated sites
 while read line; do modified=${line::-3}; echo $modified ; done < all2yag_aftermaf_ahmm.prune.in|tr ":" "\t" > pruned_sites_forvcf






# Generate sample list for bcftools query
grep "AFP\|CYP\|TX" puma_yagasref.bamlist |awk '{print $0}' > sample_list4ahmm


# prepare input vcf
bcftools convert filtered_all2yag.vcf -o filtered_all2yag.vcf.gz
bcftools index filtered_all2yag.vcf.gz


# Run bcftools query to get partial samples (TX/CYP/AFP) for A_HMM and format it
#Nine fields before samples, include header
SAMPLE_LIST=sample_list4ahmm
REGION_FILE=pruned_sites_forvcf
INPUT_VCF=filtered_all2yag.vcf.gz
OUTPUT_NAME=all2yag_pruned4ahmm.vcf
OUTPUT_NAME2=all2yag_pruned4ahmm_unformat.vcf


bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%ID[\t%GT:%AD]\n' -R $REGION_FILE -S $SAMPLE_LIST $INPUT_VCF -o $OUTPUT_NAME &




bcftools view -o $OUTPUT_NAME2 -R $REGION_FILE -S $SAMPLE_LIST $INPUT_VCF --threads 15 &

lin.yuan@tilden:/space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm$ grep "#CHROM" all2yag_pruned4ahmm_unformat.vcf
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AFP10.yag.bam   AFP11.yag.bam   AFP12.yag.bam   AFP13.yag.bam   AFP14.yag.bam   AFP15.yag.bam   AFP16.yag.bam   AFP17.yag.bam    AFP18.yag.bam   AFP19.yag.bam   AFP1.yag.bam    AFP20.yag.bam   AFP21.yag.bam   AFP22.yag.bam   AFP23.yag.bam   AFP24.yag.bam   AFP25.yag.bam   AFP26.yag.bam   AFP27.yag.bam    AFP28.yag.bam   AFP29.yag.bam   AFP2.yag.bam    AFP30.yag.bam   AFP31.yag.bam   AFP3.yag.bam    AFP4.yag.bam    AFP5.yag.bam    AFP6.yag.bam    AFP7.yag.bam    AFP8.yag.bam    AFP9.yag.bam     CYP45.yag.bam   CYP47.yag.bam   CYP51.yag.bam   CYP60.yag.bam   TX101.yag.bam   TX105.yag.bam   TX106.yag.bam   TX107.yag.bam   TX108.yag.bam










# The header of “all2yag_pruned4ahmm.vcf” has [int] and we don’t want that. We will manually modify it (delete “[integer]”).


grep "#CHROM" all2yag_pruned4ahmm_unformat.vcf > newheader
# manually delete all the spaces and rename col9


# then concatenate two things
tail -n+2 all2yag_pruned4ahmm.vcf > noheader.vcf 
cat newheader noheader.vcf > all2yag_pruned_reheaded.vcf &




# create “sample2population” file
# follow the format
lin.yuan@tilden:/space/s1/lin.yuan/puma/analysis_yagAsReference/ahmm$ more sample2population
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
AFP29.yag.bam   0
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
