#!/bin/bash
# arguments: num genes, num factors, num cells, iterations
SCODE=../../../../tools/SCODE/SCODE.jl
DATA_PATH=../../../../data/mESC/

k_vals=(2 4 6 8)
for k in ${k_vals[@]}; do
    mkdir SCODE_D_$k
    julia --project=. $SCODE $DATA_PATH/exp.txt $DATA_PATH/time.txt SCODE_D_$k 100 $k 456 1000
done
