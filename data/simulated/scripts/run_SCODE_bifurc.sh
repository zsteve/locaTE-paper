#!/bin/bash
# arguments: num genes, num factors, num cells, iterations
SCODE_DIR=/home/stephenz/stephenz/locaTE-paper/tools/SCODE

k_vals=(2 4 6 8)
N_reps=5
for i in $(seq 1 $N_reps); do 
for k in ${k_vals[@]}; do
    dir="SCODE_D_"$k
    mkdir $dir
    for j in $(ls exp_branch*.txt | awk -F'branch' '{ print $2 }' | cut -d'.' -f 1); do 
		julia --project=$SCODE_DIR $SCODE_DIR/SCODE.jl exp_branch$j.txt time_branch$j.txt $dir $(wc -l exp_branch$j.txt | cut -d' ' -f 1) $k $(wc -l time_branch$j.txt | cut -d' ' -f 1) 5
		mv $dir/A.txt $dir"/A_branch"$j"_rep_"$i".txt"
	done
done
done
