#!/bin/bash
#SBATCH --job-name="locaTE"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=mig
#SBATCH --time=0-02:00:00
#SBATCH --array=1-1
#SBATCH --output=%x-%A_%a.out

srcpath="/home/stephenz/stephenz/locaTE-paper/scripts"
datapath=/data/gpfs/projects/punim0638/stephenz/locaTE-paper/data/simulated/Synthetic_1000/dyn-BFStrange/dyn-BFStrange-1000-9

parameters=`sed -n "${SLURM_ARRAY_TASK_ID} p" $datapath/params`
parameterArray=($parameters)

k=${parameterArray[0]}
lamda1=${parameterArray[1]}
lamda2=${parameterArray[2]}

suffix=$k"_"$lamda1"_"$lamda2

JULIA=julia
mkdir $datapath/locate_output"_"$suffix
for i in $(ls -d "$datapath/P_"*".npy"); do
	ptype=$(echo $i | awk -F'P_' '{ print $2 }' | cut -d'.' -f 1)
	echo Transition matrix $ptype
	JULIA_NUM_THREADS=4 $JULIA $srcpath/infer_locate.jl --X $datapath/X.npy --X_pca $datapath/X_pca.npy --P $i --C $datapath/C.npy --k $k --lambda1 $lamda1 --lambda2 $lamda2 --outdir $datapath/locate_output"_"$suffix/ --suffix $ptype
done

