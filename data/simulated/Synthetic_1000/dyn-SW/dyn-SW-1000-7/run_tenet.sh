kvals=$(seq 1 8)
for k in ${kvals[@]}; do 
	bash ./TENET X.csv 1 dpt.csv branch.csv $k
	sleep 5
	mv TE_result_matrix.txt A_tenet_k_$k.txt
done
