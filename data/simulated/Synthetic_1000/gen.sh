#!/bin/bash
DIR_CUR=$(pwd)

for i in $(ls); do 
	for j in $(ls $i | grep -E "1000-([0-9]|10)$"); do
		DIR="$(pwd)/$i/$j"
		echo "Preprocessing $DIR"
		# python ~/stephenz/sc-causal-grn/data_benchmarking/preprocess.py $DIR --nbranches 0
		# cp run*.sh $DIR/
		# cp params* $DIR/
		# sed -i "s~__DATAPATH__~$DIR~g" $DIR/run.sh
		# sed -i "s~__DATAPATH__~$DIR~g" $DIR/run_cespgrn.sh
		# sed -i "s~__DATAPATH__~$DIR~g" $DIR/run_undir.sh
		# sbatch $DIR/run.sh
		# sbatch $DIR/run_cespgrn.sh
		# sbatch $DIR/run_undir.sh
		# Cleanup
		# cd $DIR
		# rm -r *_output_*
		# cd $DIR_CUR
		bash cleanup.sh $DIR 
	done
done
