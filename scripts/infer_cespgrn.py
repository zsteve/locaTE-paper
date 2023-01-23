import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--cespgrnpath", type = str)
parser.add_argument("--X", type = str)
parser.add_argument("--X_pca", type = str)
parser.add_argument("--k", type = int, help = "neighbourhood size")
parser.add_argument("--bandwidth", type = float, help = "bandwidth")
parser.add_argument("--lamda", type = float, help = "lambda")
parser.add_argument("--outdir", type = str, default = "./")
args = parser.parse_args()

import sys, os
sys.path.append(args.cespgrnpath)
import numpy as np
import g_admm as cg
import kernel
import torch
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

counts = np.load(args.X)
X_pca= np.load(args.X_pca)

K, K_trun = kernel.calc_kernel_neigh(X_pca, bandwidth = args.bandwidth, truncate = True, truncate_param = args.k)
empir_cov = cg.est_cov(X = counts, K_trun = K_trun, weighted_kt = True)
cespgrn = cg.G_admm_minibatch(X=counts[:, None, :], K=K, pre_cov=empir_cov, batchsize = 100)
thetas = cespgrn.train(max_iters=2500, n_intervals=100, lamb=args.lamda)

np.save(args.outdir + "G_cespgrn.npy", thetas.reshape(thetas.shape[0], -1))



