import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import os
import glob

try:
	adata = ad.read_h5ad("anndata_dpt_groups_mapped.h5ad")
except:
	adata = ad.read_h5ad("anndata.h5ad")
	sc.tl.dpt(adata)
dpt=adata.obs.dpt_pseudotime
adata = adata[np.argsort(dpt), :]
dpt = np.sort(dpt)

is_bifurc = ("dpt_groups" in adata.obs.columns)

DIR=os.getcwd()
SCRIPT_DIR="/data/gpfs/projects/punim0638/stephenz/locaTE-paper/data/simulated/Synthetic_1000/scripts/"
TOOL_DIR="/data/gpfs/projects/punim0638/stephenz/locaTE-paper/tools/"
if is_bifurc:
	TENET_SCRIPT="%s/run_tenet_bifurc.sh" % SCRIPT_DIR
	SCODE_SCRIPT="%s/run_SCODE_bifurc.sh" % SCRIPT_DIR
else:
	TENET_SCRIPT="%s/run_tenet.sh" % SCRIPT_DIR
	SCODE_SCRIPT="%s/run_SCODE.sh" % SCRIPT_DIR

# write for TENET
DIR_COMPARISON = os.path.join(DIR, "tenet/")
try:
	os.mkdir(DIR_COMPARISON)
	os.system("ln -s %s/TENET/* %s" % (TOOL_DIR, DIR_COMPARISON))
except:
	pass
pd.DataFrame(adata.X).to_csv(DIR_COMPARISON + "X.csv")
pd.DataFrame(dpt).to_csv(DIR_COMPARISON + "dpt.csv", header = False, index = False)
if is_bifurc:
	pd.DataFrame(1*np.isin(adata.obs.dpt_groups, ['0'])).to_csv(DIR_COMPARISON + "branch0.csv", header = False, index = False)
	pd.DataFrame(1*np.isin(adata.obs.dpt_groups, ['1'])).to_csv(DIR_COMPARISON + "branch1.csv", header = False, index = False)
	pd.DataFrame(1*np.isin(adata.obs.dpt_groups, ['2'])).to_csv(DIR_COMPARISON + "branch2.csv", header = False, index = False)
else:
	pd.DataFrame(np.ones(adata.shape[0], dtype = int)).to_csv(DIR_COMPARISON + "branch.csv", header = False, index = False)	

# run TENET
cmd="cd %s && ml load java && bash %s" % (DIR_COMPARISON, TENET_SCRIPT)
os.system(cmd)

# combine TENET
if is_bifurc:
	l_0 = np.isin(adata.obs.dpt_groups, ['0']).mean()
	l_1 = np.isin(adata.obs.dpt_groups, ['1']).mean()
	l_2 = np.isin(adata.obs.dpt_groups, ['2']).mean()
	for k in range(1, 8+1):
		A_0 = pd.read_csv(DIR_COMPARISON + "A_tenet_branch0_k_%d.txt" % k, sep = "\t", index_col = 0)
		A_1 = pd.read_csv(DIR_COMPARISON + "A_tenet_branch1_k_%d.txt" % k, sep = "\t", index_col = 0)
		A_2 = pd.read_csv(DIR_COMPARISON + "A_tenet_branch2_k_%d.txt" % k, sep = "\t", index_col = 0)
		(l_0*A_0 + l_1*A_1 + l_2*A_2).to_csv(DIR_COMPARISON + "A_tenet_combined_k_%d.txt" % k)

# write for SCODE
DIR_COMPARISON = os.path.join(DIR, "scode/")
try:
	os.mkdir(DIR_COMPARISON)
except:
	pass

if is_bifurc:
	branch_mask = np.isin(adata.obs.dpt_groups, ['0'])
	pd.DataFrame(adata.X[branch_mask, :].T).to_csv(DIR_COMPARISON + "exp_branch0.txt", header = False, index = False, sep = "\t")
	pd.DataFrame(dpt[branch_mask]).to_csv(DIR_COMPARISON + "time_branch0.txt", header = False, sep = "\t")
	branch_mask = np.isin(adata.obs.dpt_groups, ['1'])
	pd.DataFrame(adata.X[branch_mask, :].T).to_csv(DIR_COMPARISON + "exp_branch1.txt", header = False, index = False, sep = "\t")
	pd.DataFrame(dpt[branch_mask]).to_csv(DIR_COMPARISON + "time_branch1.txt", header = False, sep = "\t")
	branch_mask = np.isin(adata.obs.dpt_groups, ['2'])
	pd.DataFrame(adata.X[branch_mask, :].T).to_csv(DIR_COMPARISON + "exp_branch2.txt", header = False, index = False, sep = "\t")
	pd.DataFrame(dpt[branch_mask]).to_csv(DIR_COMPARISON + "time_branch2.txt", header = False, sep = "\t")
else:
	pd.DataFrame(adata.X.T).to_csv(DIR_COMPARISON + "exp.txt", header = False, index = False, sep = "\t")
	pd.DataFrame(dpt).to_csv(DIR_COMPARISON + "time.txt", header = False, sep = "\t")

# run SCODE
cmd="cd %s && bash %s" % (DIR_COMPARISON, SCODE_SCRIPT)
os.system(cmd)

# combine SCODE
if is_bifurc:
	l_0 = np.isin(adata.obs.dpt_groups, ['0']).mean()
	l_1 = np.isin(adata.obs.dpt_groups, ['1']).mean()
	l_2 = np.isin(adata.obs.dpt_groups, ['2']).mean()
	for fname in glob.glob(DIR_COMPARISON + "SCODE_D*"):
		for i in range(1, 5+1):
			A_0 = pd.read_csv(fname + "/A_branch0_rep_%d.txt" % i, sep = "\t", index_col = False, header = None)
			A_1 = pd.read_csv(fname + "/A_branch1_rep_%d.txt" % i, sep = "\t", index_col = False, header = None)
			A_2 = pd.read_csv(fname + "/A_branch2_rep_%d.txt" % i, sep = "\t", index_col = False, header = None)
			(l_0*A_0 + l_1*A_1 + l_2*A_2).to_csv(fname + "/A_combined_rep_%d.txt" % i, sep = "\t", index = None, header = None)

# Scribe
import Scribe as scr
DIR_COMPARISON = os.path.join(DIR, "scribe/")
try:
	os.mkdir(DIR_COMPARISON)
except:
	pass

model = scr.read_export.load_anndata(adata)
if is_bifurc == False:
	model.expression = pd.DataFrame(adata.X).T
	model.expression_raw = model.expression
	# model.expression.set_index(pd.Index(["g" + str(i) for i in range(1, adata.X.shape[1]+1)]), inplace = True)
	model.expression.set_index(adata.var.index, inplace = True)
model.rdi(delays=[5, 20, 40], number_of_processes=1, uniformization=False, differential_mode=False);
np.save(DIR_COMPARISON + "/G_scribe.npy", model.rdi_results["MAX"].to_numpy())
