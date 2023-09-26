import numpy as np
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('in_matrix', type = str)
parser.add_argument('out_matrix', type = str)
parser.add_argument('--noise_level', type = float, default = 0)
parser.add_argument('--srand', type = int, default = 0)
args = parser.parse_args()

np.random.seed(args.srand)

X = pd.read_csv(args.in_matrix, index_col = 0)

scale = np.sqrt(np.mean(X.to_numpy()**2))
E = scale * np.random.randn(*X.shape)
X_noisy = (1-args.noise_level)*X + args.noise_level*E

print("norm(X_noisy - X): ", np.linalg.norm((X - X_noisy).to_numpy()))
X_noisy.to_csv(args.out_matrix)
