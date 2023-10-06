#!/bin/bash
#SBATCH --job-name="static"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=mig
#SBATCH --time=0-00:30:00
#SBATCH --array=1-100
#SBATCH --output=%x-%A_%a.out

path=$(pwd)
datapath=`sed -n "${SLURM_ARRAY_TASK_ID} p" paths`

JULIA="/home/stephenz/julia-1.8.4/bin/julia"

# run julia first due to segfault issue
cd $datapath
$JULIA $path/run_pidc.jl --project=pidc_env

source ~/.bashrc && source ~/.profile
# run rest 
ml load MATLAB
conda activate py39
python $path/run_static.py

