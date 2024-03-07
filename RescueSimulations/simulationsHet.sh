#Calculate heterozygosity
for file in *vcf.gz; 
do 
python3 ./SlidingWindowHet.py $file scaffolds 1000000 100000 1 &
done

#Subset 10 individuals of each:
#Before bottleneck Np2=3000
cut -f1-2,3004-30014 everglades_additive_seed1.p2.55400.vcf.gz_het_1000kbWin_100Kbstep.txt > p2.55400.N3000.subset.het.txt

#After bottleneck Np1=5000, Np2=10
cut -f1-2,14- everglades_additive_seed1.p2.55411.vcf.gz_het_1000kbWin_100Kbstep.txt > p2.55411.N10.subset.het.txt
cut -f1-2,5004-50014 everglades_additive_seed1.p1.55411.vcf.gz_het_1000kbWin_100Kbstep.txt > p1.55411.N5000.subset.het.txt

#After rescue Np2=200 , rescue 55411
#5 years
cut -f1-2,204-2014 everglades_additive_seed1.p2.55416.vcf.gz_het_1000kbWin_100Kbstep.txt > p2.55416.N200.subset.het.txt
#10 years
cut -f1-2,204-2014 everglades_additive_seed1.p2.55421.vcf.gz_het_1000kbWin_100Kbstep.txt > p2.55421.N200.subset.het.txt
#50 years
cut -f1-2,204-2014 everglades_additive_seed1.p2.55461.vcf.gz_het_1000kbWin_100Kbstep.txt > p2.55461.N200.subset.het.txt
#100 years
cut -f1-2,204-2014 everglades_additive_seed1.p2.55511.vcf.gz_het_1000kbWin_100Kbstep.txt > p2.55511.N200.subset.het.txt
#1000
cut -f1-2,204-2014 everglades_additive_seed1.p2.56411.vcf.gz_het_1000kbWin_100Kbstep.txt > p2.56411.N200.subset.het.txt
#2000
cut -f1-2,204-2014 everglades_additive_seed1.p2.57410.vcf.gz_het_1000kbWin_100Kbstep.txt > p2.57410.N200.subset.het.txt
