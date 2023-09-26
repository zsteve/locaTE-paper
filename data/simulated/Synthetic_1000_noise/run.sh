#!/bin/bash

CURDIR=$(pwd)
noise_levels=(0 0.1 0.25 0.5 0.75)
dropout_levels=(0 0.1 0.25 0.5 0.75)

for d in $(ls | grep "dyn-*"); do
    cd $d
    for SRC_DIR in $(ls | grep "_dyn-"); do
        for s in ${noise_levels[@]}; do 
            for t in ${dropout_levels[@]}; do 
                TARG_DIR=${SRC_DIR:1}"_noise_"$s"_dropout_"$t
                echo $TARG_DIR
                cp -r $SRC_DIR $TARG_DIR
                python $CURDIR/apply_noise.py $SRC_DIR/VelocityData.csv $TARG_DIR/VelocityData.csv --noise_level $s 
                python $CURDIR/apply_dropout.py $SRC_DIR/ExpressionData.csv $TARG_DIR/VelocityData.csv $TARG_DIR/ExpressionData.csv $TARG_DIR/VelocityData.csv  --noise_level $t
            done
        done
    done
    cd $CURDIR
done
