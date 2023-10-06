import numpy as np
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('expr_matrix', type = str)
parser.add_argument('velo_matrix', type = str)
parser.add_argument('out_expr_matrix', type = str)
parser.add_argument('out_velo_matrix', type = str)
parser.add_argument('--noise_level', type = float, default = 0)
parser.add_argument('--srand', type = int, default = 0)
args = parser.parse_args()

np.random.seed(args.srand)

X = pd.read_csv(args.expr_matrix, index_col = 0)
V = pd.read_csv(args.velo_matrix, index_col = 0)

X_noisy = X.copy()
V_noisy = V.copy()
for i in range(X_noisy.shape[0]):
    q = np.quantile(X.iloc[i, :], args.noise_level)
    idx=(X_noisy.iloc[i, :] < q) & (np.random.rand(X_noisy.shape[1]) < args.noise_level)
    X_noisy.iloc[i, idx] = 0 
    V_noisy.iloc[i, idx] = 0 

print("norm(X_noisy - X): ", np.linalg.norm((X - X_noisy).to_numpy()))
print("norm(V_noisy - V): ", np.linalg.norm((V - V_noisy).to_numpy()))

X_noisy.to_csv(args.out_expr_matrix)
V_noisy.to_csv(args.out_velo_matrix)
