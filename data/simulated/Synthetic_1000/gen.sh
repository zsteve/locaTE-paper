#!/bin/bash
DIR_CUR=$(pwd)

declare -A branch
branch=( ["dyn-BF"]=1 ["dyn-BFC"]=0 ["dyn-BFStrange"]=1 ["dyn-CY"]=0 ["dyn-LI"]=0 ["dyn-LL"]=0 ["dyn-SW"]=0 ["dyn-TF"]=1 )

for i in $(ls); do 
	for j in $(ls $i | grep -E "1000-([0-9]|10)$"); do
		DIR="$(pwd)/$i/$j"
		echo "Preprocessing $DIR"
		python ../../../scripts/preprocess_boolode.py $DIR --nbranches ${branch[$i]}
		# cp scripts/run*.sh $DIR/ # copy all run scripts 
		# cp scripts/params* $DIR/ # copy param sets
		# sed -i "s~__DATAPATH__~$DIR~g" $DIR/run.sh
		# sed -i "s~__DATAPATH__~$DIR~g" $DIR/run_cespgrn.sh
		# sed -i "s~__DATAPATH__~$DIR~g" $DIR/run_undir.sh
		# echo "Submitting batch job"
		# sbatch $DIR/run.sh
		# sbatch $DIR/run_cespgrn.sh
		# sbatch $DIR/run_undir.sh
		# Cleanup
		# cd $DIR
		# rm -r *_output_*
		# cd $DIR_CUR
		# bash scripts/cleanup.sh $DIR 
	done
done
