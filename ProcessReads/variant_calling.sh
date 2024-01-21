ref=/space/s1/lin.yuan/puma/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
bcftools mpileup -b puma_yagRef.bamlist -a AD,DP,SP -Ou --fasta-ref $ref --threads 15 --ff 1804 --rf 2 -Q 20 -q 10 | bcftools call -a GQ -c --ploidy 2 -O z -o ./unfiltered_all2yag.vcf.gz --threads 15
