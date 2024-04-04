#!/bin/bash

for i in {1..10}
do   
  #Create directory for "chromosome"
  mkdir seed$i
  mv *seed$seed.* seed$seed
  cd seed$seed
  
  #Change name of chromosome
  echo -e "1\t$i" > chromosome_mapping_file.txt
  for file in *vcf
  do
    out=${file%.vcf.gz}_chr$i.vcf.gz
    bcftools annotate --rename-chrs chromosome_mapping_file.txt $file -Oz -o $out
    #Calculate heterozygosity
    python3 ./SlidingWindowHet.py $out scaffolds 1000000 100000 $i &
  done
  
done

