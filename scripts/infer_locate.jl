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

mi_all = zeros(size(X, 1), size(X, 2)^2);
p = Progress(size(X, 1))

P = npzread(args["P"])
P_sp = sparse(P^args["k"])
π_unif = fill(1/size(P, 1), size(P, 1))'
Q = (P' .* π_unif)./(π_unif * P)';
QT_sp = sparse((Q^args["k"])')
# directed inference
@info "Directed inference"
@info "Computing RDI scores"
alg = DiscretizeBayesianBlocks()
disc = locaTE.discretizations_bulk(X; alg = alg)

for i = 1:size(X, 1)
	mi_all[i, :] = get_MI(X, compute_coupling(X, i, P_sp, QT_sp, R_sp), gene_idxs[:, 1], gene_idxs[:, 2]; disc = disc, alg = alg)
	next!(p)
end

w = vec(maximum(mi_all; dims = 2))
w ./= mean(w)
@info "Applying CLR"
mi_all_clr = apply_wclr(mi_all, size(X, 2))
mi_all_clr[isnan.(mi_all_clr)] .= 0
@info "Denoising"
G = fitsp(mi_all_clr, L; λ1 = args["lambda1"], λ2 = args["lambda2"], maxiter = Nq)

npzwrite(string(args["outdir"], "G_$(args["suffix"]).npy"), G)
npzwrite(string(args["outdir"], "mi_$(args["suffix"]).npy"), mi_all)
npzwrite(string(args["outdir"], "mi_clr_$(args["suffix"]).npy"), mi_all_clr)
