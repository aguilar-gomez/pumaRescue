ls /space/s2/diana/puma/allbams_Pyaguarundi/bam_output_allsampletoOutgroup/CYP*bam > CFPbamfiles 
ls /space/s2/diana/puma/allbams_Pyaguarundi/bam_output_allsampletoOutgroup/AFP*bam > PTFPbamfiles 
ls /space/s2/diana/puma/allbams_Pyaguarundi/bam_output_allsampletoOutgroup/TX*bam > TXbamfiles 


#Get sites 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' puma_simplePASS_variants_all.vcf.gz  > puma_variants
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' puma_simplePASS_all.vcf.gz  > puma_sites
SITES=puma_sites
angsd index sites $SITES
#Create Saf files
#No SNP filtering, keep invariable
bamlist=$1
REF=
angsd -bam $bamlist -out $SITES -minMapQ 25 -minQ 25 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -GL 1 -doMaf 1 -doMajorMinor 4 -doGlf 2 -nThreads 10 -ref $REF -setMaxDepthInd 50 -doCounts 1 -skipTriallelic 1

#SNP filtering
#Do all filters and generate beagle simultaneously
/home/diana/bin/angsd -bam $bamlist -out pumilio_version11 -minInd 250 -setMinDepthInd 1 \ 
                      -minMapQ 25 -minQ 25 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
                      -GL 1 -doMaf 1 -doMajorMinor 4 -doGlf 2 -SNP_pval 1e-6 -nThreads 20 \
                      -ref $REF -setMaxDepthInd 50 -doCounts 1 -skipTriallelic 1 -dosnpstat 1 \
                      -doHWE 1 -sb_pval 1e-4 -hetbias_pval 1e-6 -doGeno 3 -sites pumilio_sites_all.pos \
                      -rf pumilio_sites_all.rf -doBcf 1 -doPost 1 --ignore-RG 0 -edge_pval 1e-4 \
                      -mapQ_pval 1e-4





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
