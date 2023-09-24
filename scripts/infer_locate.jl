using locaTE
using OptimalTransport
using NPZ
using StatsBase
using SparseArrays
using ProgressMeter
using NearestNeighbors
using Graphs
using GraphSignals
using Printf
using ArgParse
using LinearAlgebra
using Discretizers
using NNlib

s = ArgParseSettings()
@add_arg_table s begin
    "--X"
        help = "path to counts matrix"
        arg_type = String
    "--X_pca"
        help = "path to PCA coordinates"
        arg_type = String
    "--P"
        help = "path to transition matrix"
        arg_type = String
    "--C"
        help = "path to cost matrix"
        arg_type = String
    "--eps"
        help = "ε for QOT neighbourhood construction"
        arg_type = Float64
        default = 2.5
    "--k"
        help = "power for transition matrix"
        arg_type = Int
        default = 1
    "--k_lap"
        help = "number of neighbours for Laplacian"
        arg_type = Int 
        default = 15
    "--lambda1"
        arg_type = Float64
        default = 5.0
    "--lambda2"
        arg_type = Float64
        default = 0.01
    "--outdir"
        arg_type = String
        default = "./"
    "--suffix"
        arg_type = String
        default = ""
end

args = parse_args(s)
Nq = 500

@info "Loading files"
X = npzread(args["X"])
# for Bayesian blocks, cutoff artifactually small counts 
X = relu.(X .- 10^(-0.5))
X_pca = npzread(args["X_pca"])
C = npzread(args["C"])

@info "Constructing neighbourhood kernel" 
R = quadreg(ones(size(X, 1)), ones(size(X, 1)), C, args["eps"]*mean(C))
R_sp = sparse(R)

gene_idxs = vcat([[j, i]' for i = 1:size(X, 2) for j = 1:size(X, 2)]...);

# construct kNN and Laplacian
kdtree = KDTree(X_pca')
idxs, dists = knn(kdtree, X_pca', args["k_lap"]);
A = spzeros(size(X_pca, 1), size(X_pca, 1));
for (i, j) in enumerate(idxs)
    A[i, j] .= 1.0
end
L = sparse(normalized_laplacian(max.(A, A'), Float64));

P = npzread(args["P"])
P_sp = sparse(P^args["k"])
Q = to_backward_kernel(P);
QT_sp = sparse((Q^args["k"])')
# directed inference
@info "Directed inference"
@info "Computing TE scores"
alg = DiscretizeBayesianBlocks()
disc = locaTE.discretizations_bulk(X; alg = alg)

TE = estimate_TE(X, 1:size(X, 2), 1:size(X, 2), P_sp, QT_sp, R_sp)
w = vec(maximum(TE; dims = 2))
w ./= mean(w)
@info "Applying CLR"
TE_clr = apply_wclr(TE, size(X, 2))
TE_clr[isnan.(TE_clr)] .= 0
# @info "Denoising"
G = fitsp(TE_clr, L; λ1 = args["lambda1"], λ2 = args["lambda2"], maxiter = Nq)

# CDF normalization
A = reshape(maximum(G; dims = 1), size(X, 2), size(X, 2))
G_cdf = apply_cdf_norm(G, A .+ 1e-9);
# Static
cdf_norm2(x) = cdf_norm(x, x .+ 1e-9)
G_static=reshape(mean(G; dims = 1), size(X, 2), size(X, 2))
G_static_cdf=cdf_norm2(G_static)

npzwrite(string(args["outdir"], "G_$(args["suffix"]).npy"), G)
npzwrite(string(args["outdir"], "G_cdf_$(args["suffix"]).npy"), G_cdf)
npzwrite(string(args["outdir"], "G_static_$(args["suffix"]).npy"), G_static)
npzwrite(string(args["outdir"], "G_static_cdf_$(args["suffix"]).npy"), G_static_cdf)
npzwrite(string(args["outdir"], "TE_$(args["suffix"]).npy"), TE)
npzwrite(string(args["outdir"], "TE_clr_$(args["suffix"]).npy"), TE_clr)
