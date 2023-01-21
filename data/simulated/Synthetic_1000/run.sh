#!/bin/bash
#SBATCH --job-name="scgrn"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=mig
#SBATCH --time=0-03:00:00
#SBATCH --array=1-125
#SBATCH --output=%x-%A_%a.out

parameters=`sed -n "${SLURM_ARRAY_TASK_ID} p" params`
parameterArray=($parameters)

k=${parameterArray[0]}
lamda1=${parameterArray[1]}
lamda2=${parameterArray[2]}

suffix=$k"_"$lamda1"_"$lamda2
srcpath="/data/gpfs/projects/punim0638/stephenz/sc-causal-grn/src"
datapath=__DATAPATH__

JULIA=/home/stephenz/julia-1.8.0/bin/julia
# ml load julia 
mkdir $datapath/infer_output"_"$suffix
for i in $(ls -d "$datapath/P_"*".npy"); do
# for i in $(ls -d "$datapath/P_"*".npy" | grep -e dpt -e statot -e pba); do # only run on dpt, statot, pba kernels 
# for i in $(ls -d "$datapath/P_"*".npy" | grep -e pba); do # only PBA 
# for i in $(ls -d "$datapath/P_"*".npy" | grep -e statot_ent); do # only statot_ent
	ptype=$(echo $i | awk -F'P_' '{ print $2 }' | cut -d'.' -f 1)
	echo Transition matrix $ptype
	$JULIA $srcpath/infer.jl --X $datapath/X.npy --X_pca $datapath/X_pca.npy --P $i --C $datapath/C.npy --k $k --lambda1 $lamda1 --lambda2 $lamda2 --outdir $datapath/infer_output"_"$suffix/ --suffix $ptype
done

