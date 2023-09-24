using Pkg
Pkg.activate("sincerities_env")
Pkg.add(["DataFrames", "NPZ", "StatsBase", "LinearAlgebra", "NNlib", "XLSX", "Discretizers"])
using DataFrames
using NPZ
using StatsBase
using LinearAlgebra
using NNlib
using XLSX
using Discretizers

idx = Colon()
X = npzread("X.npy")[:, idx]
dpt = npzread("dpt.npy")
try
    global genes = CSV.read("genes.txt", Array)[:, 2][idx];
catch e
    global genes = ["gene_$i" for i = 1:size(X, 2)];
end

ld = LinearDiscretizer(binedges(DiscretizeUniformWidth(10), dpt))
dpt_bins = encode(ld, dpt)
XLSX.writetable("X.xlsx", insertcols!(DataFrame(X, genes), size(X, 2)+1, :t => dpt_bins), overwrite = true)
