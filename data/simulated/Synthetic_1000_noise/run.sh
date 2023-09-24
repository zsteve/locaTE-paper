#!/bin/bash

CURDIR=$(pwd)

for d in $(ls | grep "dyn-*"); do
    cd $d
    SRC_DIR=$(ls | grep "_dyn-")
    noise_levels=(0 0.1 0.25 0.5 1.0)
    for s in ${noise_levels[@]}; do 
        TARG_DIR=${SRC_DIR:1}"_noise_"$s
        echo $TARG_DIR
        cp -r $SRC_DIR $TARG_DIR
        python $CURDIR/apply_noise.py $SRC_DIR/VelocityData.csv $TARG_DIR/VelocityData.csv --level $s
    done
    cd $CURDIR
done
