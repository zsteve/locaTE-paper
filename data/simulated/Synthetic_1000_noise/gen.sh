#!/bin/bash
DIR_CUR=$(pwd)

for i in $(ls); do 
    for j in $(ls $i | grep -E "_noise_"); do
        DIR="$(pwd)/$i/$j"
        echo "Preprocessing $DIR"
        # python ../../../scripts/preprocess_boolode.py $DIR 
        cp scripts/run*.sh $DIR/ # copy all run scripts 
        cp scripts/params* $DIR/ # copy param sets
        sed -i "s~__DATAPATH__~$DIR~g" $DIR/run.sh
        # sed -i "s~__DATAPATH__~$DIR~g" $DIR/run_cespgrn.sh
        # sed -i "s~__DATAPATH__~$DIR~g" $DIR/run_undir.sh
        echo "Submitting batch job"
        sbatch $DIR/run.sh
        # sbatch $DIR/run_cespgrn.sh
        # sbatch $DIR/run_undir.sh
    done
done
