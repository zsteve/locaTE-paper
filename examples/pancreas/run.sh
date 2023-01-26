#!/bin/bash
DATA_DIR="../../data/pancreas/"
julia ../../tools/locaTE.jl/src/locaTE_cmd.jl --lambda1 10.0 --lambda2 0.001 --k_lap 25 --suffix locate --gpu $DATA_DIR/X.npy $DATA_DIR/X_pca.npy $DATA_DIR/P.npy $DATA_DIR/R.npy
