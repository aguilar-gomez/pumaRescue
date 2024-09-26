#!/bin/bash
#Job name:
#SBATCH --job-name=hetero
#
# Account:
#SBATCH --account=co_moorjani
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
#SBATCH --time=240:00:00
#
##SBATCH --mem=300G
## Command(s) to run :
for i in {1..10}
do
	/global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i -d m_rate=0.27 evg_central_amr_GR.slim  > evg_central_amr_GR_${i}.out &
	/global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i evg_central_amr_noGR.slim  > evg_central_amr_noGR_${i}.out &
#	/global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i -d m_rate=0.27 cfp_tx_GR.slim  > cfp_tx_GR_${i}.out &
#	/global/home/users/zhangyulin9806/SLiM4/build/slim -d input_seed=$i cfp_tx_noGR.slim  > cfp_tx_noGR_${i}.out &
done
wait

