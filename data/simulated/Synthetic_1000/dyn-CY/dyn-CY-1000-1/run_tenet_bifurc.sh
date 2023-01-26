kvals=$(seq 1 8)
for k in ${kvals[@]}; do 
	bash ./TENET X.csv 1 dpt.csv branch0.csv $k
	mv TE_result_matrix.txt A_tenet_branch0_k_$k.txt
	bash ./TENET X.csv 1 dpt.csv branch1.csv $k
	mv TE_result_matrix.txt A_tenet_branch1_k_$k.txt
	bash ./TENET X.csv 1 dpt.csv branch2.csv $k
	mv TE_result_matrix.txt A_tenet_branch2_k_$k.txt
done
