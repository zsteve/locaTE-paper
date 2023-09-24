import numpy as np
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('in_matrix', type = str)
parser.add_argument('out_matrix', type = str)
parser.add_argument('--level', type = float, default = 0)
args = parser.parse_args()

X = pd.read_csv(args.in_matrix, index_col = 0)
scale = np.sqrt(np.mean(X.to_numpy()**2))
E = args.level * scale * np.random.randn(*X.shape)
X_noisy = X + E
print("norm(X_noisy - X): ", np.linalg.norm((X - X_noisy).to_numpy()))
X_noisy.to_csv(args.out_matrix)
