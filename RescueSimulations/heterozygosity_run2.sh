for ((i=1; i<=10; i++)); do     echo -e "$i\t100000000"; done > scaffolds


#!/bin/bash
SCRIPT=/space/s2/diana/puma/simulation/simulation_run2/SlidingWindowHet.py
for i in {1..10}
do   
  #Create directory for "chromosome"
  mkdir -p seed$i
  mv *seed$i.* seed$i
  cp scaffolds seed$i
  cd seed$i
  
  #Change name of chromosome
  echo -e "1\t$i" > chromosome_mapping_file.txt
  for file in *vcf
  do
    bgzip $file
    tabix $file.gz
    out=${file%.vcf}_chr$i.vcf.gz
    bcftools annotate --rename-chrs chromosome_mapping_file.txt $file.gz -Oz -o $out
    #Calculate heterozygosity
    python3 $SCRIPT $out scaffolds 100000 10000 $i &
  done
  cd ..
done



