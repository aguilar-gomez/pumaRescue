
#Generate beagle file with angsd
SITES=./sites_angsd
SCAFFOLDS=./region_angsd
bamlist=./bamlist
REF=../../genome_outgroup/GCF_014898765.1_PumYag_genomic.fna

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

#3%
/home/diana/programs/ohana/tools/sample-sites.py ./puma.lgm 3 ./puma_3percent.lgm



#Do OHANA as usual using the .lgm files

# run qpas
Kmin=3 # min value of k
Kmax=6 # maximum value of k
# don’t run everything in one time; it is too thread-consuming
mi=450 # maximum number of iterations
e=0.08 # the minimum between the likelihood difference per iteration
qpas=/home/diana/programs/ohana/bin/qpas
lgmName=./puma_3percent.lgm
matrixName=puma_3pct

for k in $(seq $Kmin $Kmax)
do
echo "running qpas with k $k"
$qpas $lgmName -k $k -qo ${matrixName}_k${k}_e${e}_mi${mi}_q.matrix -fo ${matrixName}_k${k}_e${e}_mi${mi}_f.matrix -e ${e} -mi $mi > out.qpas_${matrixName}_k${k}mi${mi}e${e} &
done



# run nemeco tree
Kmin=3 # min value of k
Kmax=6 # maximum value of k
OHANA=/home/diana/programs/ohana/bin
mi=450 # maximum number of iterations
e=0.08 # the minimum between the likelihood difference per iteration
lgmName=./puma_3percent.lgm
matrixName=puma_3pct

for k in $(seq $Kmin $Kmax);
do
echo "running nemeco with k $k"
$OHANA/nemeco $lgmName ${matrixName}_k${k}_e${e}_mi${mi}_f.matrix -co c.matrix_${matrixName}_k${k} -mi 5 > out_${matrixName}_.nemk${k}
$OHANA/convert cov2nwk c.matrix_${matrixName}_k${k} ${matrixName}_k${k}.nwk
tail -n +2 ${matrixName}_k${k}_e${e}_mi${mi}_q.matrix  > ${matrixName}_k${k}.Q
done
