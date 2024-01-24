# recode to the compatible format
INPUT=all2yag_after_maf
OUTPUT_NAME=all2yag_aftermaf_ohana
plink_1.90_6.26 --bfile $INPUT --recode 12 tab --out $OUTPUT_NAME --threads 16 --allow-extra-chr




# down sampling 3%
INPUT=all2yag_aftermaf_ohana.ped
OUTPUT_DGM=all2yag.dgm
DOWN_SAMPLE_PERCENT=3
OUTPUT_NAME=all2yag_3percent.dgm
/home/diana/programs/ohana/bin/convert ped2dgm $INPUT $OUTPUT_DGM
/home/diana/programs/ohana/tools/sample-sites.py $OUTPUT_DGM $DOWN_SAMPLE_PERCENT $OUTPUT_NAME





# run qpas
Kmin=3 # min value of k
Kmax=8 # maximum value of k
# don’t run everything in one time; it is too thread-consuming
mi=450 # maximum number of iterations
e=0.08 # the minimum between the likelihood difference per iteration; actually 1 is small enough
qpas=/home/diana/programs/ohana/bin/qpas
dgmName=all2yag_3percent.dgm
matrixName=all2yag_3pct
# run ohana structure (qpas)
for k in $(seq $Kmin $Kmax)
do
echo "running qpas with k $k"
$qpas $dgmName -k $k -qo ${matrixName}_k${k}_e${e}_mi${mi}_q.matrix -fo ${matrixName}_k${k}_e${e}_mi${mi}_f.matrix -e ${e} -mi $mi > out.qpas_${matrixName}_k${k}mi${mi}e${e} &
done




# run nemeco tree
Kmin=3 # min value of k
Kmax=8 # maximum value of k
OHANA=/home/diana/programs/ohana/bin
mi=450 # maximum number of iterations
e=0.08 # the minimum between the likelihood difference per iteration; actually 1 is small enough
dgmName=all2yag_3percent.dgm
matrixName=all2yag_3pct


for k in $(seq $Kmin $Kmax);
do
echo "running nemeco with k $k"
$OHANA/nemeco $dgmName ${matrixName}_k${k}_e${e}_mi${mi}_f.matrix -co c.matrix_${matrixName}_k${k} -mi 5 > out_${matrixName}_.nemk${k}
$OHANA/convert cov2nwk c.matrix_${matrixName}_k${k} ${matrixName}_k${k}.nwk
tail -n +2 ${matrixName}_k${k}_e${e}_mi${mi}_q.matrix  > ${matrixName}_k${k}.Q
done
