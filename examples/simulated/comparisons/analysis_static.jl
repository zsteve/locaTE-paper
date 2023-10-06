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

dataset_dirs = glob([DATA_DIR, r"dyn.*-*-[0-9]*$", ])
sim = split(split(dataset_dirs[1], r"-1$")[1], "/")[end]
datasets = tmap(process_dataset, dataset_dirs)
Ng = size(first(datasets)["X"], 2)
Nc = size(first(datasets)["X"], 1)

@info "Loading job_data_static"
job_data_static = tmap(x -> load_job_data(x, Ng; locate_what = "G_static_cdf", cespgrn_what = "G_static", static = true), dataset_dirs);

is_simple = !any(occursin.(["dyn-BFStrange", "dyn-SW"], sim))

is_bifurcating = any(occursin.(["dyn-BFStrange-", "dyn-BF-", "dyn-TF-"], sim))
if is_bifurcating 
	outfiles_tenet = map(x -> glob(joinpath(x, "tenet/A_tenet_combined_k_*.txt")), dataset_dirs)
	outfiles_scode = map(x -> glob(joinpath(x, "scode/SCODE_D_4/A_combined_rep_*.txt")), dataset_dirs)
else
	outfiles_tenet = map(x -> glob(joinpath(x, "tenet/A_tenet_k_*.txt")), dataset_dirs)
	outfiles_scode = map(x -> glob(joinpath(x, "scode/SCODE_D_4/A_rep_*.txt")), dataset_dirs)
end
outfiles_scribe = map(x -> glob(joinpath(x, "scribe/G_scribe.npy")), dataset_dirs)
outfiles_pidc = map(x -> glob(joinpath(x, "G_pidc.npy")), dataset_dirs)
outfiles_sincerities = map(x -> glob(joinpath(x, "sincerities/A.txt")), dataset_dirs)
outfiles_grisli = map(x -> glob(joinpath(x, "grisli/A_grisli_L_*.csv")), dataset_dirs)

J_static = map(x -> mean(x["J"]; dims = 1)[1, :, :] .> 0.05, datasets);
baseline_static = map(x -> mean(x), J_static);

scores_static = Dict("tenet" => [map(x -> process_tenet(x, J_static[i])[2] / baseline_static[i], outfiles_tenet[i]) for i = 1:length(outfiles_tenet)],
    "scode" => [map(x -> process_scode(x, J_static[i])[2] / baseline_static[i], outfiles_scode[i]) for i = 1:length(outfiles_scode)],
    "scribe" => [map(x -> process_scribe(x, J_static[i])[2] / baseline_static[i], outfiles_scribe[i]) for i = 1:length(outfiles_scribe)], 
    "pidc" => [map(x -> process_pidc(x, J_static[i])[2] / baseline_static[i], outfiles_pidc[i]) for i = 1:length(outfiles_pidc)], 
    "sincerities" => [map(x -> process_sincerities(x, J_static[i])[2] / baseline_static[i], outfiles_sincerities[i]) for i = 1:length(outfiles_sincerities)], 
    "grisli" => [map(x -> process_grisli(x, J_static[i])[2] / baseline_static[i], outfiles_grisli[i]) for i = 1:length(outfiles_grisli)])

for k in (is_simple ? filter(x -> !any(occursin.(["symm", "statot", "pba"], x)), keys(first(job_data_static)["outputs"])) : filter(x -> !any(occursin.(["symm", "statot_ent", ], x)), keys(first(job_data_static)["outputs"])))
    scores_static[k] = [map(x -> process_locate(x, J_static[i])[2], job_data_static[i]["outputs"][k]) / baseline_static[i] for i = 1:length(job_data_static)]
end

k = "locaTE_velo_dot"
i = 9
map(x -> process_locate(x, J_static[i])[2], job_data_static[i]["outputs"][k][1:10])

plt=plot(; legend = nothing, ylabel = "AUPR ratio", title = string("Directed static networks: ", sim_name[sim]), rotation = 45, bottom_margin = 5Plots.mm, size = (PLT_CELL, PLT_CELL))
for k in filter(x -> occursin("locaTE_", x), sort(collect(keys(scores_static)); rev = true))
    boxplot!([map_label(k), ], map(NaNMath.maximum, scores_static[k]), color = :lightgrey, linecolor = :black)
end
boxplot!(["CeSpGRN", ], map(NaNMath.maximum, scores_static["cespgrn"]), color = :lightgrey, linecolor = :black)
boxplot!(["TENET", ], map(NaNMath.maximum, scores_static["tenet"]), color = :lightgrey, linecolor = :black)
boxplot!(["PIDC", ], vcat(scores_static["pidc"]...), color = :lightgrey, linecolor = :black)
boxplot!(["Scribe", ], vcat(scores_static["scribe"]...), color = :lightgrey, linecolor = :black)
boxplot!(["SCODE", ], vcat(scores_static["scode"]...), color = :lightgrey, linecolor = :black)
boxplot!(["SINCERITIES", ], vcat(scores_static["sincerities"]...), color = :lightgrey, linecolor = :black)
boxplot!(["GRISLI", ], map(NaNMath.maximum, scores_static["grisli"]), color = :lightgrey, linecolor = :black)
savefig(string(FIG_DIR, sim, "_static_aupr.pdf"))
plt

# plot representative outputs for static networks
make_square(x) = reshape(x, Int(sqrt(length(x))), Int(sqrt(length(x))))
k = 1
plt_true = heatmap(J_static[k]; c = :Greys, colorbar = nothing, title = "True");
# locaTE
i = argmax(scores_static["locaTE_velo_dot"][k])
A = make_square(mean(job_data[k]["outputs"]["locaTE_velo_dot"][i]; dims = 1))
plt_locate = heatmap(A; clim = (0, quantile(vec(A), 0.95)), c = :Greys, colorbar = nothing, title = "locaTE-velo");
# locaTE-dpt
i = argmax(scores_static["locaTE_dpt"][k])
A = make_square(mean(job_data[k]["outputs"]["locaTE_dpt"][i]; dims = 1))
plt_locate_dpt = heatmap(A; clim = (0, quantile(vec(A), 0.95)), c = :Greys, colorbar = nothing, title = "locaTE-dpt");
# CeSpGRN
i = argmax(scores_static["cespgrn"][k])
A = make_square(mean(job_data[k]["outputs"]["cespgrn"][i]; dims = 1))
plt_cespgrn = heatmap(A; clim = (0, quantile(vec(A), 0.95)), c = :Greys, colorbar = nothing, title = "CeSpGRN");
# TENET
i = argmax(scores_static["tenet"][k])
A = process_tenet(outfiles_tenet[k][i], J_static[k])[1]
plt_tenet = heatmap(A; clim = (0, quantile(vec(A), 0.95)), c = :Greys, colorbar = nothing, title = "TENET");
# PIDC
i = argmax(scores_static["pidc"][k])
A = process_pidc(outfiles_pidc[k][i], J_static[k])[1]
plt_pidc = heatmap(A; clim = (0, quantile(vec(A), 0.95)), c = :Greys, colorbar = nothing, title = "PIDC");
# Scribe
i = argmax(scores_static["scribe"][k])
A = process_scribe(outfiles_scribe[k][i], J_static[k])[1]
plt_scribe = heatmap(A; clim = (0, quantile(vec(A), 0.95)), c = :Greys, colorbar = nothing, title = "Scribe");
# SCODE
i = argmax(scores_static["scode"][k])
A = process_scode(outfiles_scode[k][i], J_static[k])[1]
plt_scode = heatmap(A; clim = (0, quantile(vec(A), 0.95)), c = :Greys, colorbar = nothing, title = "SCODE");
# SINCERITIES
i = argmax(scores_static["sincerities"][k])
A = process_sincerities(outfiles_sincerities[k][i], J_static[k])[1]
plt_sincerities = heatmap(A; clim = (0, quantile(vec(A), 0.95)), c = :Greys, colorbar = nothing, title = "SINCERITIES");
# GRISLI
i = argmax(scores_static["grisli"][k])
A = process_grisli(outfiles_grisli[k][i], J_static[k])[1]
plt_grisli = heatmap(A; c = :Greys, colorbar = nothing, title = "GRISLI");
# 
plt=plot(plt_locate, plt_locate_dpt, plt_cespgrn, plt_tenet, 
		plt_pidc, plt_scribe, plt_scode, plt_sincerities, plt_grisli, plt_true, layout = (2, 5), size = (1_000, 400))
savefig(string(FIG_DIR, sim, "_static_best.pdf"))
plt
