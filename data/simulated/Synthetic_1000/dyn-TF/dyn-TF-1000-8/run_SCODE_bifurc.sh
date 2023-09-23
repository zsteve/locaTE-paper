#!/bin/bash
# arguments: num genes, num factors, num cells, iterations
SCODE_DIR=/home/stephenz/stephenz/locaTE-paper/tools/SCODE

k_vals=(2 4 6 8)
N_reps=5
for i in $(seq 1 $N_reps); do 
for k in ${k_vals[@]}; do
    dir="SCODE_D_"$k
    mkdir $dir
    julia --project=$SCODE_DIR $SCODE_DIR/SCODE.jl exp_branch0.txt time_branch0.txt $dir $(wc -l exp_branch0.txt | cut -d' ' -f 1) $k $(wc -l time_branch0.txt | cut -d' ' -f 1) 1000
    mv $dir/A.txt $dir/A_branch0_rep_$i.txt
    julia --project=$SCODE_DIR $SCODE_DIR/SCODE.jl exp_branch1.txt time_branch1.txt $dir $(wc -l exp_branch1.txt | cut -d' ' -f 1) $k $(wc -l time_branch1.txt | cut -d' ' -f 1) 1000
    mv $dir/A.txt $dir/A_branch1_rep_$i.txt
    julia --project=$SCODE_DIR $SCODE_DIR/SCODE.jl exp_branch2.txt time_branch2.txt $dir $(wc -l exp_branch2.txt | cut -d' ' -f 1) $k $(wc -l time_branch2.txt | cut -d' ' -f 1) 1000
    mv $dir/A.txt $dir/A_branch2_rep_$i.txt
done
done
