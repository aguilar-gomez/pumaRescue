#!/bin/bash
#Job name:
#SBATCH --job-name=recgr5
#
# Account:
#SBATCH --account=fc_moorjani
#
# Partition:
#SBATCH --partition=savio3_bigmem
#
#QoS:
##SBATCH --qos=savio_normal
#
# Request one node:
#SBATCH --nodes=1
#
# Specify one task:
#SBATCH --ntasks=1
#
# Number of processors for single task needed for use case:
#SBATCH --cpus-per-task=32
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
##SBATCH --mem=390G
#SBATCH --exclude=n0006.savio3,n0007.savio3,n0008.savio3,n0009.savio3
## Command(s) to run :

nrep=5
#for i in {16..20}
#do
#	/global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i -d output_seed=$((i + nrep*20)) -d rep=$nrep cfp_tx_h05_ori_GR.slim  > ori${nrep}/cfp_tx_h05_GR_${i}.out &
#	/global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i -d output_seed=$((i + nrep*20)) -d rep=$nrep cfp_tx_ori_GR.slim  > ori${nrep}/cfp_tx_GR_${i}.out &
#done

nrep=5
#for i in {1..15}
#do
#	/global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i -d output_seed=$((i + nrep*20)) -d rep=$nrep cfp_tx_h05_ori_noGR.slim  > ori${nrep}/cfp_tx_h05_noGR_${i}.out &
#	/global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i -d output_seed=$((i + nrep*20)) -d rep=$nrep cfp_tx_ori_noGR.slim  > ori${nrep}/cfp_tx_noGR_${i}.out &
#done

nrep=5
#for i in {16..20}
#do
#	/global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i -d output_seed=$((i + nrep*20)) -d rep=$nrep cfp_tx_recentsplit_GR.slim  > recentsplit${nrep}/cfp_tx_recentsplit_GR_${i}.out &
#	/global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i -d output_seed=$((i + nrep*20)) -d rep=$nrep cfp_tx_oldsplit_GR.slim  > oldsplit${nrep}/cfp_tx_oldsplit_GR_${i}.out &
#done

nrep=5
for i in {1..15}
do
#        /global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i -d output_seed=$((i + nrep*20)) -d rep=$nrep cfp_tx_recentsplit_noGR.slim  > recentsplit${nrep}/cfp_tx_recentsplit_noGR_${i}.out &
       /global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i -d output_seed=$((i + nrep*20)) -d rep=$nrep cfp_tx_oldsplit_noGR.slim  > oldsplit${nrep}/cfp_tx_oldsplit_noGR_${i}.out &
done

wait
