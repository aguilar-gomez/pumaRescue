# prepare force component file

# 'I' Individuals and 'K' Components
50 6
# Population Assignments per Individual
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 1 1 2 2 2 2 3 3 4
4 4 4 4 5 5 5 5 5 6




#Population 0 AFP
12 1
0 0 0 0 0 0
1 1 1 1 1 1


# Population 1 BR
12 1
1 0 0 0 0 0
1 0 0 0 0 0


# Population 2 CYP
12 1
0 1 0 0 0 0
0 1 0 0 0 0


# Population 3 EVG
12 1
0 0 1 0 0 0
0 0 1 0 0 0


# Population 4 CAL (SC SMM)
12 1
0 0 0 1 0 0
0 0 0 1 0 0


# Population 5 TX
12 1
0 0 0 0 1 0
0 0 0 0 1 0


# Population 6 YNP
12 1
0 0 0 0 0 1
0 0 0 0 0 1





# copy the site files from last time’s run
cp -s /space/s1/lin.yuan/puma/analysis_yagAsReference/down_sampling/all2yag_3percent.dgm .



# run supervised qpas
k=6
mi=450 # maximum number of iterations
e=0.08 # the minimum between the likelihood difference per iteration; actually 1 is small enough
qpas=/home/diana/programs/ohana/bin/qpas
dgmName=all2yag_3percent.dgm
matrixName=supervised_3pct
force=fg.txt
$qpas $dgmName -k $k -qo ${matrixName}_k${k}_e${e}_mi${mi}_q.matrix -fo ${matrixName}_k${k}_e${e}_mi${mi}_f.matrix -e ${e} -mi $mi -fg $force 






# run nemeco tree
k=6
OHANA=/home/diana/programs/ohana/bin
mi=450 # maximum number of iterations
e=0.08 # the minimum between the likelihood difference per iteration; actually 1 is small enough
dgmName=all2yag_3percent.dgm
matrixName=supervised_3pct


$OHANA/nemeco $dgmName ${matrixName}_k${k}_e${e}_mi${mi}_f.matrix -co c.matrix_${matrixName}_k${k} -mi $mi > out_${matrixName}_.nemk${k}
$OHANA/convert cov2nwk c.matrix_${matrixName}_k${k} ${matrixName}_k${k}.nwk
tail -n+2 ${matrixName}_k${k}_e${e}_mi${mi}_q.matrix  > ${matrixName}_k${k}.Q


