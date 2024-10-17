ls /space/s2/diana/puma/allbams_Pyaguarundi/bam_output_allsampletoOutgroup/CYP*bam > CFPbamfiles 
ls /space/s2/diana/puma/allbams_Pyaguarundi/bam_output_allsampletoOutgroup/AFP*bam > PTFPbamfiles 
ls /space/s2/diana/puma/allbams_Pyaguarundi/bam_output_allsampletoOutgroup/TX*bam > TXbamfiles 


#Get sites 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' puma_simplePASS_variants_all.vcf.gz  > puma_variants
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' puma_simplePASS_all.vcf.gz  > puma_sites
cut -f1-2 puma_sites > puma_sites.txt
angsd sites index puma_sites.txt
cut -f1 puma_sites.txt|sort|uniq > puma_sites.rf

#Create Saf files
#No SNP filtering, keep invariable
#!/bin/bash
POP=$1
SITES=puma_sites.txt
RF=puma_sites.rf
bamlist=${POP}bamfiles
REF=/space/s1/lin.yuan/puma/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna
angsd -bam $bamlist -out $POP -sites $SITES -rf $RF -minMapQ 30 -minQ 25 -remove_bads 1 \
      -uniqueOnly 1 -only_proper_pairs 1 -GL 1 -nThreads 10 -ref $REF -doSaf 1 -baq 1

#Run for populations
nohup sh doSaf.sh CFP > out.CFP &

#Fst index and scan use unfolded saf
pop1=$1
pop2=$2
#joint SFS
realSFS $pop1.r.saf.idx $pop2.r.saf.idx -P 5 >$pop1.$pop2.ml
#Fst index
realSFS fst index $pop1.r.saf.idx $pop2.r.saf.idx -sfs $pop1.$pop2.ml  -fstout $pop1.$pop2
realSFS fst stats $pop1.$pop2.fst.idx > fstindex$pop1$pop2
#Fst scan
realSFS fst stats2 $pop1.$pop2.fst.idx -win 100000 -step 20000 -type 2 > fst_w100kb$pop1$pop2


#PBS code and Windows of SNPs
#PBS and fst all together
pop1=$1
pop2=$2
pop3=$3
SNP=$4
realSFS fst index $pop1.r.saf.idx $pop2.r.saf.idx $pop3.r.saf.idx -sfs $pop1.$pop2.ml -sfs $pop1.$pop3.ml -sfs $pop2.$pop3.ml -fstout $pop1.$pop2.$pop3
realSFS fst stats  $pop1.$pop2.$pop3.fst.idx>pbs.$pop1.$pop2.$pop3
realSFS fst stats2 $pop1.$pop2.$pop3.fst.idx -win 100000 -step 20000 -type 2 >slidingwindow100Ks20k.$pop1.$pop2.$pop3
#Make SNP windows
realSFS fst print $pop1.$pop2.$pop3.fst.idx  > $pop1.$pop2.$pop3.printed.txt
PBS_SNP_windows.py $pop1.$pop2.$pop3.printed.txt $pop1.$pop2.$pop3.wSNP$SNP $SNP

#Tajimas D and thetas 
#The values of the different thetas are negative because they are logscale!

#Calculate thetas
realSFS saf2theta $pop.saf.idx -sfs $pop.folded.sfs -outname $pop
thetaStat print $pop.thetas.idx > $pop.thetas.persite.txt

#For Dxy
nohup realSFS -P 4 $POP.r.saf.idx -sites $pop1.$pop2.shared.pos -fold 1 1> Results/$POP.sfs 2>Results/$POP.outsfs &

#Calculate Dxy
$NGSTOOLS/ngsStat -npop 2 -postfiles $pop1.shared.saf $pop2.shared.saf -nsites $NSITES -nind $pop1_nind $pop2_nind -outfile $pop1.$pop2.stats.txt


#Do windows of SNPS
Tajimas_SNPmidpoint.py $pop.thetas.persite.shared.sorted $pop.Tajimas.midpoint.SNP$sites $sites $popsize &
SNPwindows_midpoint.py $pop1.$pop2.persite.Dxy $pop1.$pop2.SNP200.Dxy.midpoint 200 Dxy &
