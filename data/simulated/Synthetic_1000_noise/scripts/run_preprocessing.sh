#!/bin/bash
#SBATCH --job-name="locaTE-preproc"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=mig
#SBATCH --time=0-01:00:00
#SBATCH --output=%x-%A_%a.out

srcpath="/home/stephenz/stephenz/locaTE-paper/scripts"
datapath=__DATAPATH__

source ~/.bashrc && source ~/.profile 
conda activate py39
python $srcpath/preprocess_boolode.py $datapath
