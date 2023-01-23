# Script for pre-processing BoolODE output
import argparse 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
import cellrank as cr
import scvelo as scv
import scipy as sp
import statot
import os
import networkx as nx
import sklearn
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

parser = argparse.ArgumentParser()
parser.add_argument('dir', type=str)
parser.add_argument('--nbranches', type=int, default = 0)
# point to PBA directory (https://github.com/AllonKleinLab/PBA) 
parser.add_argument('--pbadir', type=str, default="/data/gpfs/projects/punim0638/stephenz/locaTE-paper/tools/pba/") 
parser.add_argument('--nneighbors', type=int, default=25) 
args = parser.parse_args()

adata = anndata.AnnData(pd.read_csv(args.dir+"/ExpressionData.csv", index_col=0).T)
adata.obs["dpt"] = pd.read_csv(args.dir+"/PseudoTime.csv", index_col=0).max(1)
adata.obsm["v"] = pd.read_csv(args.dir+"/VelocityData.csv", index_col=0).sort_index().T
adata.obsm["J"] = pd.read_csv(args.dir+"/JacobianData.csv", index_col=0).sort_index().T
adata.obsm["X_raw"] = np.copy(adata.X)

# sc.pp.log1p(adata)
sc.pp.pca(adata)
try:
    # try doing kNN. sometimes this can fail if dropouts cause repeated cells (i.e. all zero)
    sc.pp.neighbors(adata, n_neighbors = args.nneighbors, method = 'gauss')
except:
    # kNN failed
    print("Warning: kNN failed. Injecting small amount of noise and trying again")
    adata.X += 1e-6*np.random.randn(*adata.X.shape)
    sc.pp.neighbors(adata, n_neighbors = args.nneighbors, method = 'gauss')

adata.uns['iroot'] = np.argmin(adata.obs.dpt)
sc.tl.diffmap(adata)
sc.tl.dpt(adata, n_branchings = args.nbranches)

# check if tSNE already computed
try:
    adata.obsm["X_tsne"] = pd.read_csv(args.dir+"/tsne100.tsv", index_col=0, sep = "\t")
except:
    sc.tl.tsne(adata, perplexity=100)
    plt.figure() 
    plt.scatter(adata.obsm["X_tsne"][:, 0], adata.obsm["X_tsne"][:, 1], c = adata.obs.dpt)
    plt.savefig(args.dir+"/tsne100.png")

def compute_velocity_kernel(k, sigma):
    # compute velocity-induced transition matrix with kernel k and bandwidth sigma
    if sigma is None:
        # median heuristic, see e.g. https://github.com/theislab/cellrank/blob/f6d281e62a5cee976d25973b5a5e4c1050f37429/cellrank/kernels/_velocity_kernel.py#L230 
        logits_all = np.concatenate([k(adata.obsm["v"].iloc[i, :].to_numpy().reshape(1, -1).astype(np.float64),
                        (adata.X[adata.uns['neighbors']['connectivities'][i, :].indices, :] - adata.X[i, :]).astype(np.float64),
                                softmax_scale =  1.0)[1] for i in range(adata.shape[0])])
        sigma = 1/np.median(logits_all)
    probs_all_fwd = [k(adata.obsm["v"].iloc[i, :].to_numpy().reshape(1, -1).astype(np.float64),
                    (adata.X[adata.uns['neighbors']['connectivities'][i, :].indices, :] - adata.X[i, :]).astype(np.float64),
                    softmax_scale = sigma)[0] for i in range(adata.shape[0])]
    P_velo = sp.sparse.csr_matrix(adata.uns['neighbors']['connectivities'].shape)
    for i, p in enumerate(probs_all_fwd):
        P_velo[i, adata.uns['neighbors']['connectivities'][i, :].indices] = p
    return P_velo

# apply velocity kernels 
for (k, sigma, kerneltype) in zip([cr.tl.kernels.DotProductScheme(), cr.tl.kernels.CorrelationScheme(), cr.tl.kernels.CosineScheme()],
                            [None, None, None], ["dot", "corr", "cos"]):
    adata.obsm["P_velo_%s" % kerneltype] = compute_velocity_kernel(k, sigma)
# diffusion pseudotime kernel
k = cr.tl.kernels.PseudotimeKernel(adata, time_key = "dpt_pseudotime").compute_transition_matrix(threshold_scheme = "soft") 
adata.obsm["P_dpt"] = k.transition_matrix

# construct cost matrix along discrete manifold 
### kNN
G_sp = adata.obsp["distances"]
# prune edges to prevent "short-circuit" artifacts
# NB. a better solution would be to use e.g. the discrete heat kernel (https://arxiv.org/abs/2211.00805) 
def prune_edges(G, q = 0.95):
    bc = nx.edge_current_flow_betweenness_centrality(G)
    bc_id, bc_weights = zip(*bc.items()); bc_weights = np.array(bc_weights); bc_id = list(bc_id)
    bad_edges = np.where(bc_weights > np.quantile(bc_weights, q))[0]
    bad_edge_centralities = bc_weights[bad_edges]
    # now sort
    G_pruned = G.copy()
    for (i, j) in [bc_id[x] for x in bad_edges[np.argsort(bad_edge_centralities)[::-1]]]:
        G_pruned.remove_edge(i, j)
        if not nx.is_connected(G_pruned):
            G_pruned.add_edge(i, j)
    return G_pruned, bad_edges

G = nx.from_numpy_matrix(np.array(G_sp.todense()))
G_pruned, _ = prune_edges(G)
G_sp_pruned = nx.adj_matrix(G_pruned)
adata.obsm["C"] = sp.sparse.csgraph.floyd_warshall(G_sp_pruned, directed = False)**2

# Stationary OT
sink_idx = (adata.obs.dpt >= np.quantile(adata.obs.dpt, 0.925)) & (adata.obs.dpt < np.quantile(adata.obs.dpt, 0.975))
R = np.zeros(adata.shape[0])
R[sink_idx] = -50/sum(sink_idx)
R[~sink_idx] = -R.sum()/sum(~sink_idx)
gamma, mu, nu = statot.inference.statot(adata.obsm["X_pca"], g = np.exp(R), dt = 1, C = adata.obsm["C"], eps = 0.05*adata.obsm["C"].mean(), method = "quad")
adata.obsm["P_statot"] = statot.inference.row_normalise(gamma)
# entropic
gamma, mu, nu = statot.inference.statot(adata.obsm["X_pca"], g = np.exp(R), dt = 1, C = adata.obsm["C"], eps = 0.001*adata.obsm["C"].mean(), method = "ent")
gamma[gamma < np.quantile(gamma, 0.95, axis=-1).reshape(-1, 1)] = 0 
adata.obsm["P_statot_ent"] = statot.inference.row_normalise(gamma)

# PBA
os.system("cd %s && mkdir -p pba && cd pba && ln -s %s/* ." % (args.dir, args.pbadir))
np.save(args.dir + "/pba/X.npy", adata.X)
with open(args.dir + '/pba/edges.csv', 'w+') as f:
    G_pba = G_pruned # nx.from_numpy_matrix(np.array(adata.obsp["distances"].todense()))
    f.write('\n'.join('%d,%d' % x for x in G_pba.edges))
np.save(args.dir + "/pba/R.npy", R)
# to run PBA, need to use a python2.7 conda env
# os.system("cd %s/pba && conda activate py27 && pwd && \
#             python PBA_pipeline.py -e edges.csv -R R.npy -D 1.0 && cp P.npy ../P_pba.npy" % args.dir)
os.system("cd %s/pba && \
            conda run -n py27 python PBA_pipeline.py -e edges.csv -R R.npy -D 1.0 && cp P.npy ../P_pba.npy" % args.dir)

# plot transition matrices
try:
    x_ord = np.argsort(adata.obsm["X_tsne"].to_numpy()[:, 1])
except:
    x_ord = np.argsort(adata.obsm["X_tsne"][:, 1])
for k in [x for x in adata.obsm.keys() if "P_" in x]:
    try:
        P = np.array(adata.obsm[k].todense())
    except:
        P = np.array(adata.obsm[k])
    plt.imshow(np.linalg.matrix_power(P, 3)[x_ord, :][:, x_ord], vmax = np.quantile(P, 0.99))
    plt.title("%s" % k)
    plt.savefig(args.dir + "/P_%s.png" % k)

# save h5ad 
adata.write_h5ad(args.dir + "/anndata.h5ad")

# write numpy files
np.save(args.dir + "/X.npy", adata.X)
np.save(args.dir + "/X_raw.npy", adata.obsm["X_raw"])
np.save(args.dir + "/X_pca.npy", adata.obsm["X_pca"])
np.save(args.dir + "/X_tsne.npy", adata.obsm["X_tsne"])
np.save(args.dir + "/X_diffmap.npy", adata.obsm["X_diffmap"])

# save transition matrices
for k in [x for x in adata.obsm.keys() if "P_" in x]:
    print("Writing transition matrix %s..." % k)
    try:
        np.save(args.dir + "/%s.npy" % k, np.array(adata.obsm[k].todense()))
    except:
        np.save(args.dir + "/%s.npy" % k, np.array(adata.obsm[k]))

np.save(args.dir + "/C.npy", adata.obsm["C"])
np.save(args.dir + "/dpt.npy", adata.obs["dpt"].to_numpy())
np.save(args.dir + "/J.npy", adata.obsm["J"].to_numpy().reshape(-1, adata.X.shape[1], adata.X.shape[1]))
