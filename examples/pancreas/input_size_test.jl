using Distributed
@everywhere begin
	using Pkg
	Pkg.activate(".")
	using locaTE
	using LinearAlgebra
end
using NPZ
using SparseArrays
using StatsBase
using Random
using NNlib
using ProgressMeter

DATA_DIR = "../../data/pancreas/"
X = npzread(string(DATA_DIR, "X.npy"))
R = sparse(npzread(string(DATA_DIR, "R.npy")))
P = sparse(npzread(string(DATA_DIR, "P.npy")))
QT = to_backward_kernel(P);
@everywhere begin
	X=$X
	P=$P
	QT=$QT
	R=$R
end

# estimate_TE(X, 1:3, 1:3, P, QT, R)

@everywhere begin
	i = 1
	disc = locaTE.discretizations_bulk(X);
	targets, regulators = 1:size(X, 2), 1:size(X, 2)
	gene_idxs = vcat([[j, i]' for i in targets for j in regulators]...);
	clusters = I(size(P, 1))
	clusters_norm = convert(Matrix{eltype(P)}, clusters)
	clusters_norm ./= sum(clusters_norm; dims = 1);
end

TE = @showprogress pmap(i -> get_MI(X, compute_coupling(X, clusters_norm[:, i], P, QT, R), gene_idxs[:, 1], gene_idxs[:, 2]; disc = disc), 
			1:size(X, 1))
