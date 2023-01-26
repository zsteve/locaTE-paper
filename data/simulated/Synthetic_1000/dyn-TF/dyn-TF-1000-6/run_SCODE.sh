#!/bin/bash
# arguments: num genes, num factors, num cells, iterations
SCODE=/home/stephenz/stephenz/SCODE/SCODE.jl
SCODE_DIR=/home/stephenz/stephenz/SCODE

k_vals=(2 4 6 8)
for i in $(seq 1 $N_reps); do 
for k in ${k_vals[@]}; do
    dir="SCODE_D_"$k
    mkdir $dir
    julia --project=$SCODE_DIR $SCODE exp.txt time.txt $dir $(wc -l exp.txt | cut -d' ' -f 1) $k $(wc -l time.txt | cut -d' ' -f 1) 1000
    mv $dir/A.txt $dir/A_rep_$i.txt
done
done
