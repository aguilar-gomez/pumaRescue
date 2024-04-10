#!/bin/sh
#$ -l highp,h_rt=72:00:00,h_data=24G
#$ -pe shared 4
#$ -N Relatedness
#$ -cwd
#$ -m bea
#$ -o ./relate.out
#$ -e ./relate.err
#$ -M daguilar

NGSRELATE=~/bin/ngsRelate/ngsRelate
$NGSRELATE -h puma_simplePASS_variants_all.vcf.gz -O relatedness_puma -T GT


#nind:50 overall_number_of_sites:3402662

#Add names of samples
readarray -t items < "samples.list"
output_file="combinations_output.txt"
echo -e "indva\tindvb" >> "$output_file"
length=${#items[@]}
for ((i = 0; i < length; i++)); do
	for ((j = i + 1; j < length; j++)); do
    	echo -e "${items[i]}\t${items[j]}" >> "$output_file"
	done
done

paste combinations_output.txt relatedness_puma > relatedness_puma_header.txt


#Add options in nsgRelate
#!/bin/sh
#$ -l highp,h_rt=72:00:00,h_data=24G
#$ -pe shared 4
#$ -N Relatedness
#$ -cwd
#$ -m bea
#$ -o ./relate.out
#$ -e ./relate.err
#$ -M daguilar

NGSRELATE=~/bin/ngsRelate/ngsRelate
# l is minMaf

$NGSRELATE -h puma_simplePASS_variants_all.vcf.gz -O relatedness_puma_maf5 -T GT -l .05 -z samples.list -p 10 -e 0.05
