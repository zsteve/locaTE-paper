kvals=$(seq 1 16)
for k in ${kvals[@]}; do 
	bash ./TENET X.csv 16 dpt.csv cellmask.csv $k
	mv TE_result_matrix.txt A_tenet_$k.txt
done
