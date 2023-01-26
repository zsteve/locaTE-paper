#!/bin/bash
#SBATCH --job-name="cespgrn"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=mig
#SBATCH --time=0-00:30:00
#SBATCH --array=1-1
#SBATCH --output=%x-%A_%a.out

srcpath="/home/stephenz/stephenz/locaTE-paper/scripts"
datapath=/home/stephenz/stephenz/locaTE-paper/data/simulated/Synthetic_1000/dyn-BFC/dyn-BFC-1000-6

parameters=`sed -n "${SLURM_ARRAY_TASK_ID} p" $datapath/params_cespgrn`
parameterArray=($parameters)

k=${parameterArray[0]}
bandwidth=${parameterArray[1]}
lamda=${parameterArray[2]}

suffix=$k"_"$bandwidth"_"$lamda

source /home/stephenz/.bash_profile && conda activate py37
mkdir $datapath/cespgrn_output"_"$suffix
python $srcpath/cespgrn.py --cespgrnpath /home/stephenz/stephenz/locaTE-paper/misc/CeSpGRN/src --X $datapath/X.npy --X_pca $datapath/X_pca.npy --k $k --bandwidth $bandwidth --lamda $lamda --outdir $datapath/cespgrn_output"_"$suffix/
