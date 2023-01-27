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

dataset_dirs = glob([DATA_DIR, r"dyn.*-1000-[0-9]*$", ])
sim = split(split(dataset_dirs[1], r"-1$")[1], "/")[end]
datasets = tmap(process_dataset, dataset_dirs)
Ng = size(first(datasets)["X"], 2)
Nc = size(first(datasets)["X"], 1)

job_data = tmap(x -> load_job_data(x), dataset_dirs);

job_scores = tmap(((x, y), ) -> score_job_outputs(x["outputs"], y), zip(job_data, datasets));

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

sim_name = Dict("dyn-SW-1000" => "Switch", 
                "dyn-BFStrange-1000" => "Bifurcating",
                "dyn-LI-1000" => "Linear",
                "dyn-LL-1000" => "Linear long",
                "dyn-BF-1000" => "Bifurcating simple",
                "dyn-CY-1000" => "Cyclic",
                "dyn-TF-1000" => "Trifurcating",
                "dyn-BFC-1000" => "Bifurcating cycle")

is_simple = !any(occursin.(["dyn-BFStrange", "dyn-SW"], sim))
if is_simple
    # remove StatOT 
    df = df[:, filter(x -> !any(occursin.(["statot", "pba"], x)), names(df))];
else
    # remove StatOT_ent
    df = df[:, filter(x -> !any(occursin.(["statot_ent", ], x)), names(df))];
end

# best over parameter sweep
plt1=plot(title = string("Directed cell-specific networks: ", sim_name[sim]), ylabel = "AUPR ratio", legend = true, bottom_margin = 5Plots.mm, size = (PLT_CELL, PLT_CELL))
for (i, label) in enumerate(filter(x -> ~occursin("symm", x), names(df)))
    boxplot!(plt1, [map_label(label)], vec(hcat((skipmissing(df[:, label]) |> collect)...)); color = :lightgrey, linecolor = :black, rotation = 45, label = nothing)
end
hline!([mean(perf_static ./ baseline), ]; label = "Static network")
savefig(string(FIG_DIR, sim, "_dynamic_directed_aupr.pdf"))
plt1

plt2=plot(title = string("Undirected cell-specific networks: ", sim_name[sim]), ylabel = "AUPR ratio", legend = true, bottom_margin = 5Plots.mm, size = (PLT_CELL, PLT_CELL))
for (i, label) in enumerate(filter(x -> occursin("symm", x), names(df)))
    boxplot!(plt2, [map_label(label)], vec(hcat((skipmissing(df[:, label]) |> collect)...)); color = :lightgrey, linecolor = :black, rotation = 45, label = nothing)
end
hline!([mean(perf_static_symm ./ baseline_symm), ]; label = "Static network")
savefig(string(FIG_DIR, sim, "_dynamic_undirected_aupr.pdf"))
plt2

plt_12 = plot(plot(plt1; xticks = nothing), plt2; layout = (2, 1))
savefig(string(FIG_DIR, sim, "_dynamic_combined_aupr.pdf"))
plt_12

# now look at fixed parameter combinations
# idx = findfirst(x -> x == [3.0, 25.0, 0.001], job_data[1]["params"]["infer"]) # locaTE
# idx_cespgrn = findfirst(x -> x == [100.0, 0.5, 0.005], job_data[1]["params"]["cespgrn"]) # CeSpGRN
# plt2 = plot(title = "Fixed parameters", ylabel = "AUPR ratio", legend = nothing)
# boxplot!(plt2, ["locaTE"], map(x -> x["perf"][idx], job_scores)/baseline_dict["perf"]; color = :lightgrey)
# boxplot!(plt2, ["CeSpGRN"], map(x -> x["perf_cespgrn"][idx_cespgrn], job_scores)/baseline_dict["perf_cespgrn"]; color = :lightgrey)
# boxplot!(plt2, ["locaTE (symm)"], map(x -> x["perf_symm"][idx], job_scores)/baseline_dict["perf_symm"]; color = :lightgrey)

# plt=plot(plt1, plt2; size = (3/2*PLT_CELL, 2PLT_CELL/3), plot_title = "Switch, 1000 cells", left_margin = 5*Plots.mm)
# savefig(string(FIG_DIR, split(dataset_dirs[1], r"-1$")..., ".pdf"))
# plt
# plot(plt1, plt2; bottom_margin = 5Plots.mm, plot_title = string("Dynamic networks: ", split(dataset_dirs[1], r"-1$")...), titlefontsize = 12)

using CSV

is_bifurcating = any(occursin.(["dyn-BFStrange-", "dyn-BF-"], sim))
if is_bifurcating 
	outfiles_tenet = map(x -> glob(joinpath(x, "tenet/A_tenet_combined_k_*.txt")), dataset_dirs);
	outfiles_scode = map(x -> glob(joinpath(x, "scode/SCODE_D_4/A_combined_rep_*.txt")), dataset_dirs);
else
	outfiles_tenet = map(x -> glob(joinpath(x, "tenet/A_tenet_k_*.txt")), dataset_dirs);
	outfiles_scode = map(x -> glob(joinpath(x, "scode/SCODE_D_4/A_rep_*.txt")), dataset_dirs);
end
outfiles_scribe = map(x -> glob(joinpath(x, "scribe/G_scribe.npy")), dataset_dirs);
outfiles_pidc = map(x -> glob(joinpath(x, "G_pidc.npy")), dataset_dirs);

J_static = map(x -> mean(x["J"]; dims = 1)[1, :, :] .> 0.05, datasets);
baseline_static = map(x -> mean(x), J_static);

scores_static = Dict("tenet" => [map(x -> process_tenet(x, J_static[i])[2] / baseline_static[i], outfiles_tenet[i]) for i = 1:length(outfiles_tenet)],
    "scode" => [map(x -> process_scode(x, J_static[i])[2] / baseline_static[i], outfiles_scode[i]) for i = 1:length(outfiles_scode)],
    "scribe" => [map(x -> process_scribe(x, J_static[i])[2] / baseline_static[i], outfiles_scribe[i]) for i = 1:length(outfiles_scribe)], 
    "pidc" => [map(x -> process_pidc(x, J_static[i])[2] / baseline_static[i], outfiles_pidc[i]) for i = 1:length(outfiles_pidc)])
for k in (is_simple ? filter(x -> !any(occursin.(["symm", "statot", "pba"], x)), keys(first(job_data)["outputs"])) : filter(x -> !any(occursin.(["symm", "statot_ent", ], x)), keys(first(job_data)["outputs"])))
    scores_static[k] = [map(x -> process_locate(x, J_static[i])[2], job_data[i]["outputs"][k]) / baseline_static[i] for i = 1:length(job_data)]
end

plt=plot(; legend = nothing, ylabel = "AUPR ratio", title = string("Directed static networks: ", sim_name[sim]), rotation = 45, bottom_margin = 5Plots.mm, size = (PLT_CELL, PLT_CELL))
for k in filter(x -> occursin("locaTE_", x), sort(collect(keys(scores_static)); rev = true))
    boxplot!([map_label(k), ], map(NaNMath.maximum, scores_static[k]), color = :lightgrey, linecolor = :black)
end
boxplot!(["CeSpGRN", ], map(NaNMath.maximum, scores_static["cespgrn"]), color = :lightgrey, linecolor = :black)
boxplot!(["TENET", ], map(NaNMath.maximum, scores_static["tenet"]), color = :lightgrey, linecolor = :black)
boxplot!(["PIDC", ], vcat(scores_static["pidc"]...), color = :lightgrey, linecolor = :black)
boxplot!(["Scribe", ], vcat(scores_static["scribe"]...), color = :lightgrey, linecolor = :black)
boxplot!(["SCODE", ], vcat(scores_static["scode"]...), color = :lightgrey, linecolor = :black)
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
plt=plot(plt_locate, plt_locate_dpt, plt_cespgrn, plt_tenet, plt_pidc, plt_scribe, plt_scode, plt_true, layout = (2, 4), size = (800, 420))
savefig(string(FIG_DIR, sim, "_static_best.pdf"))
plt

#= uncomment to plot parameter dependences (in case of full grid search)
p = hcat(job_data[k]["params"]["locaTE_velo_dot"]...)
p_unique = unique.(eachrow(p))
y = job_scores_map[k]["locaTE_velo_dot"]

plt=plot(hline!(scatter(p_unique[1], vec(maximum(y; dims = (2, 3))), xlabel = "k", size = (3/2*PLT_CELL, PLT_CELL/2)), [perf_static, ]), 
     hline!(scatter(p_unique[2], vec(maximum(y; dims = (1, 3))), xlabel = "λ1", size = (3/2*PLT_CELL, PLT_CELL/2)), [perf_static, ]),
     hline!(scatter(p_unique[3], vec(maximum(y; dims = (1, 2))), xlabel = "λ2", size = (3/2*PLT_CELL, PLT_CELL/2)), [perf_static, ]), 
    plot_title = "locaTE", ylim = (0.9*minimum(y), 1.1*maximum(y)); legend = nothing, ylabel = "AUPRC", layout = (1, 3), bottom_margin = 5Plots.mm, left_margin = 5Plots.mm)
savefig(string(FIG_DIR, "locaTE_params_", split(dataset_dirs[1], r"-1$")..., ".pdf"))
plt

p_symm = hcat(job_data[k]["params"]["infer_velo_dot"]...)
p_symm_unique = unique.(eachrow(p_symm))
y_symm = job_scores_map[k]["infer_velo_dot_symm"]

plt=plot(hline!(scatter(p_symm_unique[1], vec(maximum(y_symm; dims = (2, 3))), xlabel = "k", size = (3/2*PLT_CELL, PLT_CELL/2)), [perf_static_symm, ]),
    hline!(scatter(p_symm_unique[2], vec(maximum(y_symm; dims = (1, 3))), xlabel = "λ1", size = (3/2*PLT_CELL, PLT_CELL/2)), [perf_static_symm, ]),
    hline!(scatter(p_symm_unique[3], vec(maximum(y_symm; dims = (1, 2))), xlabel = "λ2", size = (3/2*PLT_CELL, PLT_CELL/2)), [perf_static_symm, ]), 
    plot_title = "locaTE-symm", ylim = (0.9*minimum(y_symm), 1.1*maximum(y_symm)); legend = nothing, ylabel = "AUPRC", layout = (1, 3), bottom_margin = 5Plots.mm, left_margin = 5Plots.mm)
savefig(string(FIG_DIR, "locaTE_params_", split(dataset_dirs[1], r"-1$")..., ".pdf"))
plt

p_cespgrn = hcat(job_data[k]["params"]["cespgrn"]...)
p_cespgrn_unique = unique.(eachrow(p_cespgrn))
y_cespgrn = job_scores_map[k]["cespgrn_symm"]

plt=plot(hline!(scatter(p_cespgrn_unique[1], vec(maximum(y_cespgrn; dims = (2, 3))), xlabel = "k", size = (3/2*PLT_CELL, PLT_CELL/2)), [perf_static_symm, ]),
    hline!(scatter(p_cespgrn_unique[2], vec(maximum(y_cespgrn; dims = (1, 3))), xlabel = "bw", size = (3/2*PLT_CELL, PLT_CELL/2)), [perf_static_symm, ]),
    hline!(scatter(p_cespgrn_unique[3], vec(maximum(y_cespgrn; dims = (1, 2))), xlabel = "λ", size = (3/2*PLT_CELL, PLT_CELL/2)), [perf_static_symm, ]), 
    title = "CeSpGRN", ylim = (0.9*minimum(y_cespgrn), 1.1*maximum(y_cespgrn)); legend = nothing, ylabel = "AUPRC", layout = (1, 3), bottom_margin = 5Plots.mm, left_margin = 5Plots.mm)
savefig(string(FIG_DIR, "CeSpGRN_params_", split(dataset_dirs[1], r"-1$")..., ".pdf"))
plt

=#
