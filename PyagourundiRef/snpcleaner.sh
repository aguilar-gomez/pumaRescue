zcat unfiltered_all2yag.vcf.gz | /home/diana/bin/snpCleaner.pl -H 1e-6 -h 0 -b 1e-100 -S 1e-4 -f 1e-6 -e 1e-4 -B filtered_all2yag.bed -p all2yag_failedsites.txt.bz > filtered_all2yag.vcf
