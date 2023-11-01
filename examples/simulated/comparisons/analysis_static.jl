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
using NaNStatistics

include("process_outputs_util.jl")

s = ArgParseSettings()
@add_arg_table! s begin
	"data_dir"
		arg_type = String
		required=true
	"fig_dir"
		arg_type = String
		default = "./"
end
args = parse_args(ARGS, s)

PLT_CELL = 500
FIG_DIR=args["fig_dir"]
DATA_DIR=args["data_dir"]

# dataset_dirs = glob([DATA_DIR, r"dyn.*-*-[0-9]*$", ])
dataset_dirs = glob([DATA_DIR, r".*-*-[0-9]*$", ])
sim = split(split(dataset_dirs[1], r"-1$")[1], "/")[end]
datasets = tmap(process_dataset, dataset_dirs)
Ng = size(first(datasets)["X"], 2)
Nc = size(first(datasets)["X"], 1)

@info "Loading job_data_static"
job_data_static = tmap(x -> load_job_data(x, Ng; locate_what = "G_static_cdf", cespgrn_what = "G_static", static = true), dataset_dirs);

is_simple = !any(occursin.(["dyn-BFStrange", "dyn-SW"], sim))

function map_label(x)
    # replace method name
    l = replace(x, "_symm" => "", "locate" => "locaTE", "locate_" => "locaTE", "cespgrn" => "CeSpGRN", "velo_" => "")
end

is_bifurcating = any(occursin.(["dyn-BFStrange-", "dyn-BF-", "dyn-TF-"], sim))
if is_bifurcating 
	outfiles_tenet = map(x -> glob(joinpath(x, "tenet/A_tenet_combined_k_*.txt")), dataset_dirs)
	outfiles_scode = map(x -> glob(joinpath(x, "scode/SCODE_D_4/A_combined_rep_*.txt")), dataset_dirs)
else
	outfiles_tenet = map(x -> glob(joinpath(x, "tenet/A_tenet_k_*.txt")), dataset_dirs)
	outfiles_scode = map(x -> glob(joinpath(x, "scode/SCODE_D_4/A_rep_*.txt")), dataset_dirs)
end
outfiles_scribe = map(x -> glob(joinpath(x, "scribe/G_scribe_rdi.npy")), dataset_dirs)
outfiles_pidc = map(x -> glob(joinpath(x, "G_pidc.npy")), dataset_dirs)
outfiles_sincerities = map(x -> glob(joinpath(x, "sincerities/A.txt")), dataset_dirs)
outfiles_genie3 = map(x -> glob(joinpath(x, "genie3/A.csv")), dataset_dirs)
outfiles_grisli = map(x -> glob(joinpath(x, "grisli/A_grisli_L_*.csv")), dataset_dirs)

J_static = map(x -> mean(x["J"]; dims = 1)[1, :, :] .> 0.05, datasets);
# set all static networks diag to zero
# for x in J_static
#     x[diagind(x)] .= 0
# end
J_symm_static = map(x -> mean(x["J_symm"]; dims = 1)[1, :, :] .> 0.05, datasets);
baseline_static = mean.(J_static);
baseline_symm_static = mean.(J_symm_static); 

scores_static = Dict("tenet" => [map(x -> process_tenet(x, J_static[i])[2] / baseline_static[i], outfiles_tenet[i]) for i = 1:length(outfiles_tenet)],
    "scode" => [map(x -> process_scode(x, J_static[i])[2] / baseline_static[i], outfiles_scode[i]) for i = 1:length(outfiles_scode)],
    "scribe" => [map(x -> process_scribe(x, J_static[i])[2] / baseline_static[i], outfiles_scribe[i]) for i = 1:length(outfiles_scribe)], 
    "pidc" => [map(x -> process_pidc(x, J_static[i])[2] / baseline_static[i], outfiles_pidc[i]) for i = 1:length(outfiles_pidc)], 
    "sincerities" => [map(x -> process_sincerities(x, J_static[i])[2] / baseline_static[i], outfiles_sincerities[i]) for i = 1:length(outfiles_sincerities)], 
    "genie3" => [map(x -> process_genie3(x, J_static[i])[2] / baseline_static[i], outfiles_genie3[i]) for i = 1:length(outfiles_genie3)], 
    "grisli" => [map(x -> process_grisli(x, J_static[i])[2] / baseline_static[i], outfiles_grisli[i]) for i = 1:length(outfiles_grisli)])
for k in (is_simple ? filter(x -> !any(occursin.(["symm", "statot", "pba"], x)), keys(first(job_data_static)["outputs"])) : filter(x -> !any(occursin.(["symm", "statot_ent", ], x)), keys(first(job_data_static)["outputs"])))
    scores_static[k] = [map(x -> process_locate(x, J_static[i])[2], job_data_static[i]["outputs"][k]) / baseline_static[i] for i = 1:length(job_data_static)]
end

scores_symm_static = Dict("tenet" => [map(x -> process_tenet(x, J_symm_static[i]; symm = true)[2] / baseline_symm_static[i], outfiles_tenet[i]) for i = 1:length(outfiles_tenet)],
    "scode" => [map(x -> process_scode(x, J_symm_static[i]; symm = true)[2] / baseline_symm_static[i], outfiles_scode[i]) for i = 1:length(outfiles_scode)],
    "scribe" => [map(x -> process_scribe(x, J_symm_static[i]; symm = true)[2] / baseline_symm_static[i], outfiles_scribe[i]) for i = 1:length(outfiles_scribe)], 
    "pidc" => [map(x -> process_pidc(x, J_symm_static[i])[2] / baseline_symm_static[i], outfiles_pidc[i]) for i = 1:length(outfiles_pidc)], 
    "sincerities" => [map(x -> process_sincerities(x, J_symm_static[i]; symm = true)[2] / baseline_symm_static[i], outfiles_sincerities[i]) for i = 1:length(outfiles_sincerities)], 
    "genie3" => [map(x -> process_genie3(x, J_symm_static[i])[2] / baseline_symm_static[i], outfiles_genie3[i]) for i = 1:length(outfiles_genie3)], 
    "grisli" => [map(x -> process_grisli(x, J_symm_static[i]; symm = true)[2] / baseline_symm_static[i], outfiles_grisli[i]) for i = 1:length(outfiles_grisli)])
for k in (is_simple ? filter(x -> !any(occursin.(["symm", "statot", "pba"], x)), keys(first(job_data_static)["outputs"])) : filter(x -> !any(occursin.(["symm", "statot_ent", ], x)), keys(first(job_data_static)["outputs"])))
    scores_symm_static[k] = [map(x -> process_locate(x, J_symm_static[i]; symm = true)[2], job_data_static[i]["outputs"][k]) / baseline_symm_static[i] for i = 1:length(job_data_static)]
end

plt=plot(; legend = nothing, ylabel = "AUPR ratio", title = string("Directed static networks: ", sim_name[rsplit(sim, "-"; limit = 2)[1]]), rotation = 45, bottom_margin = 5Plots.mm, size = (PLT_CELL, PLT_CELL))
# for k in filter(x -> occursin("locaTE_", x), sort(collect(keys(scores_static)); rev = true))
for k in vcat(filter(x -> occursin("locaTE_", x) & !occursin("undir", x), sort(collect(keys(scores_static)); rev = true)), "locaTE_undir")
    boxplot!([map_label(k), ], map(NaNMath.maximum, scores_static[k]), color = :lightgrey, linecolor = :black)
end
boxplot!(["CeSpGRN", ], map(NaNMath.maximum, scores_static["cespgrn"]), color = :lightgrey, linecolor = :black)
boxplot!(["TENET", ], map(NaNMath.maximum, scores_static["tenet"]), color = :lightgrey, linecolor = :black)
boxplot!(["PIDC", ], vcat(scores_static["pidc"]...), color = :lightgrey, linecolor = :black)
boxplot!(["Scribe", ], vcat(scores_static["scribe"]...), color = :lightgrey, linecolor = :black)
boxplot!(["SCODE", ], vcat(scores_static["scode"]...), color = :lightgrey, linecolor = :black)
boxplot!(["SINCERITIES", ], vcat(scores_static["sincerities"]...), color = :lightgrey, linecolor = :black)
boxplot!(["GENIE3", ], vcat(scores_static["genie3"]...), color = :lightgrey, linecolor = :black)
boxplot!(["GRISLI", ], map(NaNMath.maximum, scores_static["grisli"]), color = :lightgrey, linecolor = :black)
savefig(string(FIG_DIR, sim, "_static_aupr.pdf"))
plt

plt=plot(; legend = nothing, ylabel = "AUPR ratio", title = string("Undirected static networks: ", sim_name[rsplit(sim, "-"; limit = 2)[1]]), rotation = 45, bottom_margin = 5Plots.mm, size = (PLT_CELL, PLT_CELL))
# for k in filter(x -> occursin("locaTE_", x), sort(collect(keys(scores_symm_static)); rev = true))
for k in vcat(filter(x -> occursin("locaTE_", x) & !occursin("undir", x), sort(collect(keys(scores_static)); rev = true)), "locaTE_undir")
    boxplot!([map_label(k), ], map(NaNMath.maximum, scores_symm_static[k]), color = :lightgrey, linecolor = :black)
end
boxplot!(["CeSpGRN", ], map(NaNMath.maximum, scores_symm_static["cespgrn"]), color = :lightgrey, linecolor = :black)
boxplot!(["TENET", ], map(NaNMath.maximum, scores_symm_static["tenet"]), color = :lightgrey, linecolor = :black)
boxplot!(["PIDC", ], vcat(scores_symm_static["pidc"]...), color = :lightgrey, linecolor = :black)
boxplot!(["Scribe", ], vcat(scores_symm_static["scribe"]...), color = :lightgrey, linecolor = :black)
boxplot!(["SCODE", ], vcat(scores_symm_static["scode"]...), color = :lightgrey, linecolor = :black)
boxplot!(["SINCERITIES", ], vcat(scores_symm_static["sincerities"]...), color = :lightgrey, linecolor = :black)
boxplot!(["GENIE3", ], vcat(scores_symm_static["genie3"]...), color = :lightgrey, linecolor = :black)
boxplot!(["GRISLI", ], map(NaNMath.maximum, scores_symm_static["grisli"]), color = :lightgrey, linecolor = :black)
savefig(string(FIG_DIR, sim, "_symm_static_aupr.pdf"))
plt

_nanargmax(x) = isnan(x[argmax(x)]) ? nanargmax(x) : argmax(x)

# plot representative outputs for static networks
k = 1
plt_true = heatmap(J_static[k]; c = :Greys, colorbar = nothing, title = "True");
# locaTE-velo-dot
i = _nanargmax(scores_static["locaTE_velo_dot"][k])
A = job_data_static[k]["outputs"]["locaTE_velo_dot"][i]
plt_locate_velo_dot = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "locaTE-velo-dot");
# locaTE-velo-cos
# i = _nanargmax(scores_static["locaTE_velo_cos"][k])
# A = job_data_static[k]["outputs"]["locaTE_velo_cos"][i]
# plt_locate_velo_cos = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "locaTE-velo-cos");
# locaTE-dpt
i = _nanargmax(scores_static["locaTE_dpt"][k])
A = job_data_static[k]["outputs"]["locaTE_dpt"][i]
plt_locate_dpt = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "locaTE-dpt");
# locaTE-undir
i = _nanargmax(scores_static["locaTE_undir"][k])
A = job_data_static[k]["outputs"]["locaTE_undir"][i]
plt_locate_undir = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "locaTE-undir");
# CeSpGRN
i = _nanargmax(scores_static["cespgrn"][k])
A = job_data_static[k]["outputs"]["cespgrn"][i]
plt_cespgrn = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "CeSpGRN");
# TENET
i = _nanargmax(scores_static["tenet"][k])
A = process_tenet(outfiles_tenet[k][i], J_static[k])[1]
plt_tenet = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "TENET");
# PIDC
i = _nanargmax(scores_static["pidc"][k])
A = process_pidc(outfiles_pidc[k][i], J_static[k])[1]
plt_pidc = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "PIDC");
# Scribe
i = _nanargmax(scores_static["scribe"][k])
A = process_scribe(outfiles_scribe[k][i], J_static[k])[1]
plt_scribe = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "Scribe");
# SCODE
i = _nanargmax(scores_static["scode"][k])
A = process_scode(outfiles_scode[k][i], J_static[k])[1]
plt_scode = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "SCODE");
# SINCERITIES
i = _nanargmax(scores_static["sincerities"][k])
A = process_sincerities(outfiles_sincerities[k][i], J_static[k])[1]
plt_sincerities = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "SINCERITIES");
# GENIE3
i = _nanargmax(scores_static["genie3"][k])
A = process_genie3(outfiles_genie3[k][i], J_static[k])[1]
plt_genie3 = heatmap(A; clim = (0, quantile(vec(A), 0.99)), c = :Greys, colorbar = nothing, title = "GENIE3");
# GRISLI
i = _nanargmax(scores_static["grisli"][k])
A = process_grisli(outfiles_grisli[k][i], J_static[k])[1]
plt_grisli = heatmap(A; c = :Greys, colorbar = nothing, title = "GRISLI");
# 
plt=plot(plt_locate_velo_dot, plt_locate_dpt, plt_locate_undir, plt_cespgrn, plt_tenet, 
		plt_pidc, plt_scribe, plt_scode, plt_sincerities, plt_genie3, plt_grisli, plt_true, layout = (3, 4), size = (1_000, 750))
savefig(string(FIG_DIR, sim, "_static_best.pdf"))
plt
