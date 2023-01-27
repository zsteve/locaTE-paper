#!/bin/bash
#SBATCH --job-name="figures"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=mig
#SBATCH --time=0-00:30:00
#SBATCH --array=1-8
#SBATCH --output=%x-%A_%a.out

path=$(pwd)
datapath=`sed -n "${SLURM_ARRAY_TASK_ID} p" paths`

JULIA="/home/stephenz/julia-1.8.4/bin/julia"
source ~/.bashrc && source ~/.bash_profile
$JULIA analysis.jl $datapath

