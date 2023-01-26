#!/bin/bash

julia ../../tools/locaTE.jl/src/locaTE_cmd.jl --lambda1 10.0 --lambda2 0.001 --k_lap 25 --suffix cmd --gpu X.npy X_pca.npy P.npy R.npy
