using Pkg
Pkg.activate(".")

import locaTE
import locaTE.prec_rec_rate 
import locaTE.aupr
using Glob
using NPZ
using OptimalTransport
using StatsBase
using Plots
using LinearAlgebra
using ProgressMeter
using NaNMath
using Suppressor
using ThreadTools
using DataFrames
using ArgParse 
using StatsPlots
using CSV

include("process_outputs_util.jl")

# s = ArgParseSettings()
# @add_arg_table! s begin
# 	"data_dir"
# 		arg_type = String
# 		required=true
# 	"fig_dir"
# 		arg_type = String
# 		default = "./"
# end
# args = parse_args(ARGS, s)

PLT_CELL = 500
# FIG_DIR=args["fig_dir"]
FIG_DIR="./"
# DATA_DIR=args["data_dir"]
DATA_DIR="../../../data/simulated/Synthetic_1000/dyn-SW/"
# DATA_DIR="../../../data/simulated/Synthetic_1000/dyn-BFStrange/"

dataset_dirs = glob([DATA_DIR, r"dyn.*-*-[0-9]*$", ])
sim = split(split(dataset_dirs[1], r"-1$")[1], "/")[end]

@info "Loading datasets"
datasets = @showprogress map(process_dataset, dataset_dirs)
Ng = size(first(datasets)["X"], 2)
Nc = size(first(datasets)["X"], 1)

@info "Loading job_data"
job_data = @showprogress map(x -> load_job_data(x, Ng; locate_what = "G_cdf", cespgrn_what = "G"), dataset_dirs);
@info "Loading job_data_static"
job_data_static = @showprogress map(x -> load_job_data(x, Ng; locate_what = "G_static_cdf", cespgrn_what = "G_static"), dataset_dirs);

@info "Scoring jobs"
job_scores = @showprogress map(((x, y), ) -> score_job_outputs(x["outputs"], y), zip(job_data, datasets));

# performance baselines
d = first(datasets)
# baseline: scGRN is same as the static true network
perf_static_symm = aupr(collect(eachcol(prec_rec_rate(d["R"]*d["J_symm_mat"], repeat(mean(d["J_symm_mat"]; dims = 1), Nc), Nq)))...)
perf_static = aupr(collect(eachcol(prec_rec_rate(d["R"]*d["J_mat"], repeat(mean(d["J_mat"]; dims = 1), Nc), Nq)))...)

job_scores_map = [score_map(x, y) for (x, y) in zip(job_scores, job_data)];
perf_best = Dict(k => map(x -> NaNMath.maximum(x[k]), job_scores) for k in keys(first(job_scores)));

# compute random predictor baselines)
compute_baseline(x, thresh) = mean(x .> thresh)
baseline = [compute_baseline(d["R"] * d["J_mat"], 0.5) for d in datasets]
baseline_symm = [compute_baseline(d["R"] * d["J_symm_mat"], 0.5) for d in datasets];

using DataFrames
using StatsPlots
baseline_dict = [Dict(k => occursin("symm", k) ? b_symm : b for k in keys(first(job_scores))) for (i, (b, b_symm)) in enumerate(zip(baseline, baseline_symm))];

df = DataFrame(perf_best) ./ vcat([DataFrame(b) for b in baseline_dict]...)
df = df[:, sort(names(df); rev = true)];

function map_label(x)
    # replace method name
    l = replace(x, "_symm" => "", "locate" => "locaTE", "locate_" => "locaTE", "cespgrn" => "CeSpGRN", "velo_" => "")
end

is_simple = !any(occursin.(["dyn-BFStrange", "dyn-SW"], sim))
if is_simple
    # remove StatOT 
    df = df[:, filter(x -> !any(occursin.(["statot", "pba"], x)), names(df))];
else
    # remove StatOT_ent
    df = df[:, filter(x -> !any(occursin.(["statot_ent", ], x)), names(df))];
end

# force undir to be the last locaTE method considered
df = df[:, [:locaTE_velo_dot_symm, :locaTE_velo_dot,
   :locaTE_velo_cos_symm, :locaTE_velo_cos,
   :locaTE_velo_corr_symm, :locaTE_velo_corr, 
   :locaTE_statot_symm, :locaTE_statot,
   :locaTE_pba_symm, :locaTE_pba,
   :locaTE_dpt_symm, :locaTE_dpt, 
   :locaTE_undir_symm, :locaTE_undir, 
   :cespgrn_symm, :cespgrn]]

# best over parameter sweep: directed cell-specific networks
plt1=plot(title = string("Directed cell-specific networks: ", sim_name[rsplit(sim, "-"; limit = 2)[1]]), ylabel = "AUPR ratio", legend = true, bottom_margin = 5Plots.mm, size = (PLT_CELL, PLT_CELL))
for (i, label) in enumerate(filter(x -> ~occursin("symm", x), names(df)))
    boxplot!(plt1, [map_label(label)], vec(hcat((skipmissing(df[:, label]) |> collect)...)); color = :lightgrey, linecolor = :black, rotation = 45, label = nothing)
end
hline!([mean(perf_static ./ baseline), ]; label = "Static network")
savefig(joinpath(FIG_DIR, string(sim, "_dynamic_directed_aupr.pdf")))
plt1

# best over parameter sweep: undirected cell-specific networks
plt2=plot(title = string("Undirected cell-specific networks: ", sim_name[rsplit(sim, "-"; limit = 2)[1]]), ylabel = "AUPR ratio", legend = true, bottom_margin = 5Plots.mm, size = (PLT_CELL, PLT_CELL))
for (i, label) in enumerate(filter(x -> occursin("symm", x), names(df)))
    boxplot!(plt2, [map_label(label)], vec(hcat((skipmissing(df[:, label]) |> collect)...)); color = :lightgrey, linecolor = :black, rotation = 45, label = nothing)
end
hline!([mean(perf_static_symm ./ baseline_symm), ]; label = "Static network")
savefig(joinpath(FIG_DIR, string(sim, "_dynamic_undirected_aupr.pdf")))
plt2

plt_12 = plot(plot(plt1; xticks = nothing), plt2; layout = (2, 1))
savefig(string(FIG_DIR, sim, "_dynamic_combined_aupr.pdf"))
plt_12
