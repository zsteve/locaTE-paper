import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import os
import glob
import sys

GENIE3_PATH="/home/stephenz/stephenz/locaTE-paper/tools/GENIE3/GENIE3_python"
sys.path.append(GENIE3_PATH)
import GENIE3

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
SCRIPT_DIR="/data/gpfs/projects/punim0638/stephenz/locaTE-paper/data/simulated/scripts/"
TOOL_DIR="/data/gpfs/projects/punim0638/stephenz/locaTE-paper/tools/"
if is_bifurc:
	TENET_SCRIPT="%s/run_tenet_bifurc.sh" % SCRIPT_DIR
	SCODE_SCRIPT="%s/run_SCODE_bifurc.sh" % SCRIPT_DIR
else:
	TENET_SCRIPT="%s/run_tenet.sh" % SCRIPT_DIR
	SCODE_SCRIPT="%s/run_SCODE.sh" % SCRIPT_DIR

# write for GENIE3 
DIR_COMPARISON = os.path.join(DIR, "genie3/")
try:
	os.mkdir(DIR_COMPARISON)
except:
	pass
pd.DataFrame(GENIE3.GENIE3(adata.X), index = adata.var.index, columns = adata.var.index).to_csv(os.path.join(DIR_COMPARISON, "A.csv"))

# write for TENET
DIR_COMPARISON = os.path.join(DIR, "tenet/")
try:
	os.mkdir(DIR_COMPARISON)
	os.system("ln -s %s/TENET/* %s" % (TOOL_DIR, DIR_COMPARISON))
except:
	pass
pd.DataFrame(adata.X).to_csv(os.path.join(DIR_COMPARISON, "X.csv"))
pd.DataFrame(dpt).to_csv(os.path.join(DIR_COMPARISON, "dpt.csv"), header = False, index = False)
if is_bifurc:
    for i in adata.obs.dpt_groups.unique():
        pd.DataFrame(1*np.isin(adata.obs.dpt_groups, [i])).to_csv(os.path.join(DIR_COMPARISON, f"branch{i}.csv"), header = False, index = False)
else:
	pd.DataFrame(np.ones(adata.shape[0], dtype = int)).to_csv(os.path.join(DIR_COMPARISON, "branch.csv"), header = False, index = False)	
# run TENET
cmd="cd %s && ml load Java && bash %s" % (DIR_COMPARISON, TENET_SCRIPT)
os.system(cmd)

# combine TENET
if is_bifurc:
	l_all = {i : np.isin(adata.obs.dpt_groups, [i]).mean()  for i in adata.obs.dpt_groups.unique()}
	for k in range(1, 8+1):
		A_all = {i : pd.read_csv(os.path.join(DIR_COMPARISON, f"A_tenet_branch{i}_k_{k}.txt"), sep = "\t", index_col = 0) for i in adata.obs.dpt_groups.unique()}
		sum([l_all[i]*A_all[i] for i in adata.obs.dpt_groups.unique()]).to_csv(os.path.join(DIR_COMPARISON, "A_tenet_combined_k_%d.txt" % k))

# write for SCODE
DIR_COMPARISON = os.path.join(DIR, "scode/")
try:
	os.mkdir(DIR_COMPARISON)
except:
	pass

if is_bifurc:
	for i in adata.obs.dpt_groups.unique():
		branch_mask = np.isin(adata.obs.dpt_groups, [i])
		pd.DataFrame(adata.X[branch_mask, :].T).to_csv(os.path.join(DIR_COMPARISON, f"exp_branch{i}.txt"), header = False, index = False, sep = "\t")
		pd.DataFrame(dpt[branch_mask]).to_csv(os.path.join(DIR_COMPARISON, f"time_branch{i}.txt"), header = False, sep = "\t")
else:
	pd.DataFrame(adata.X.T).to_csv(os.path.join(DIR_COMPARISON, "exp.txt"), header = False, index = False, sep = "\t")
	pd.DataFrame(dpt).to_csv(os.path.join(DIR_COMPARISON, "time.txt"), header = False, sep = "\t")

# run SCODE
cmd="cd %s && bash %s" % (DIR_COMPARISON, SCODE_SCRIPT)
os.system(cmd)

# combine SCODE
if is_bifurc:
	l_all = {i : np.isin(adata.obs.dpt_groups, [i]).mean()  for i in adata.obs.dpt_groups.unique()}
	for fname in glob.glob(os.path.join(DIR_COMPARISON, "SCODE_D*")):
		for j in range(1, 5+1):
			A_all = {i : pd.read_csv(os.path.join(fname, f"A_branch{i}_rep_{j}.txt"), sep = "\t", index_col = None, header = None) for i in adata.obs.dpt_groups.unique()}
			sum([l_all[i]*A_all[i] for i in adata.obs.dpt_groups.unique()]).to_csv(os.path.join(fname, "A_combined_rep_%d.txt" % j), sep = "\t", index = None, header = None)

# write for SINCERITIES 
DIR_COMPARISON = os.path.join(DIR, "sincerities/")
try:
	os.mkdir(DIR_COMPARISON)
	os.system("ln -s %s/SINCERITIES/matlab/* %s" % (TOOL_DIR, DIR_COMPARISON))
except:
	pass
_df=pd.DataFrame(adata.X, columns = adata.var.index)
_df.insert(adata.shape[1], "t", np.digitize(dpt, np.histogram_bin_edges(dpt, bins = 10), right = True))
_df.to_excel(os.path.join(DIR_COMPARISON, "X.xlsx"), header = False, index = False)
cmd="cd %s && ml load MATLAB && matlab -nodesktop -r \"run run_sincerities.m\"" % DIR_COMPARISON
os.system(cmd)

# write for GRISLI
DIR_COMPARISON = os.path.join(DIR, "grisli/")
try:
	os.mkdir(DIR_COMPARISON)
	os.system("ln -s %s/GRISLI/* %s" % (TOOL_DIR, DIR_COMPARISON))
except:
	pass
pd.DataFrame(adata.X).to_csv(os.path.join(DIR_COMPARISON, "X.csv"), header = False, index = False)
pd.DataFrame(dpt).to_csv(os.path.join(DIR_COMPARISON, "dpt.csv"), header = False, index = False)
# run GRISLI
cmd="cd %s && ml load MATLAB && matlab -nodesktop -r \"run main.m\"" % DIR_COMPARISON
os.system(cmd)

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
A_scribe_rdi = model.rdi_results["MAX"].to_numpy()
A_scribe_rdi[A_scribe_rdi == -np.inf] = 0
A_scribe = scr.Scribe.CLR(A_scribe_rdi)
np.save(os.path.join(DIR_COMPARISON, "G_scribe_rdi.npy"), A_scribe_rdi)
np.save(os.path.join(DIR_COMPARISON, "G_scribe.npy"), A_scribe)
