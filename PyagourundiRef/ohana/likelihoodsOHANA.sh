
#Generate beagle file with angsd
SITES=
SCAFFOLDS=
bamlist=
/home/diana/bin/angsd -bam $bamlist -out puma_likelihoods -minInd 50 \
        -setMinDepthInd 1 -minMapQ 25 -minQ 25 -remove_bads 1 -uniqueOnly 1 \
        -only_proper_pairs 1  -GL 1 -doMaf 1 -doMajorMinor 4 -doGlf 2 -SNP_pval 1e-4 \
        -nThreads 8 -ref $REF -setMaxDepthInd 50 -doCounts 1 -skipTriallelic 1 \
        -dosnpstat 1 -doHWE 1 -sb_pval 1e-4 -doGeno 3 -hetbias_pval 1e-6 \
        -sites $SITES -rf $SCAFFOLDS -doBcf 1 -doPost 1 \
        --ignore-RG 0 -edge_pval 1e-4 -mapQ_pval 1e-4

gunzip puma_likelihoods.beagle.gz
#OHANA
OHANA=/home/diana/programs/ohana/bin
nohup $OHANA/convert bgl2lgm ./puma_likelihoods.beagle ./puma.lgm &

#1%
/home/diana/programs/ohana/tools/sample-sites.py 

#Do OHANA as usual using the .lgm files
