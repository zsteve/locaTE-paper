kvals=$(seq 1 8)
for k in ${kvals[@]}; do 
    for i in $(ls branch*.csv | awk -F'branch' '{ print $2 }' | cut -d'.' -f 1); do 
        bash ./TENET X.csv 1 dpt.csv branch$i.csv $k
        mv TE_result_matrix.txt "A_tenet_branch"$i"_k_"$k".txt"
    done 
done
