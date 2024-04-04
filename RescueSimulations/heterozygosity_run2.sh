#!/bin/bash

for i in {1..10}
do   
  #Create directory for "chromosome"
  mkdir -p seed$i
  mv *seed$i.* seed$i
  cd seed$i
  
  #Change name of chromosome
  echo -e "1\t$i" > chromosome_mapping_file.txt
  for file in *vcf
  do
    out=${file%.vcf.gz}_chr$i.vcf
    bcftools annotate --rename-chrs chromosome_mapping_file.txt $file -Oz -o $out
    #Calculate heterozygosity
    python3 ./SlidingWindowHet.py $out scaffolds 1000000 100000 $i &
  done
  
done

