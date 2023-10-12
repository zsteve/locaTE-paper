sim_name = Dict("dyn-SW" => "Switch", 
                "dyn-BFStrange" => "Bifurcating",
                "dyn-LI" => "Linear",
                "dyn-LL" => "Linear long",
                "dyn-BF" => "Bifurcating simple",
                "dyn-CY" => "Cyclic",
                "dyn-TF" => "Trifurcating",
                "dyn-BFC" => "Bifurcating cycle")

rm_diag(A) = A - Diagonal(diag(A))

function preprocess_cespgrn(x)
    if size(x) == (Ng, Ng)
        return abs.(rm_diag(x))
    else
        return abs.(hcat([vec(rm_diag(reshape(x, Ng, Ng))) for x in eachrow(x)]...)')
    end
end

Nq = 128

function process_dataset(dir)
    J = abs.(permutedims(npzread(joinpath(dir, "J.npy")), [1, 3, 2]));
    J_symm = permutedims(cat([max.(x, x') for x in eachslice(J; dims = 1)]...; dims = 3), [3, 1, 2]);
    C = npzread(joinpath(dir, "C.npy"));
    X = npzread(joinpath(dir, "X.npy"));
    R = quadreg(ones(size(X, 1)), ones(size(X, 1)), C, 2.5*mean(C));
    dpt = npzread(joinpath(dir, "dpt.npy"))
    J_mat = reshape(J, size(J, 1), :);
    J_symm_mat = reshape(J_symm, size(J, 1), :);
    return Dict("J" => J,
                "J_symm" => J_symm,
                "C" => C,
                "X" => X,
                "R" => R,
                "dpt" => dpt,
                "J_mat" => J_mat,
                "J_symm_mat" => J_symm_mat)
end

function load_job_data(dir, Ng; locate_prefix = "locate", cespgrn_prefix = "cespgrn", locate_what = "G_cdf", cespgrn_what = "G", static = false)
    files = Dict("locaTE_velo_dot" => glob(joinpath(dir, "$(locate_prefix)_output*")), 
                "locaTE_velo_cos"  => glob(joinpath(dir, "$(locate_prefix)_output*")), 
                "locaTE_velo_corr" => glob(joinpath(dir, "$(locate_prefix)_output*")),
                "locaTE_dpt"       => glob(joinpath(dir, "$(locate_prefix)_output*")),
                "locaTE_statot"    => glob(joinpath(dir, "$(locate_prefix)_output*")),
                "locaTE_statot_ent"    => glob(joinpath(dir, "$(locate_prefix)_output*")),
                "locaTE_pba"       => glob(joinpath(dir, "$(locate_prefix)_output*")),
                "locaTE_undir"       => glob(joinpath(dir, "$(locate_prefix)_output*")),
                 "cespgrn"        => glob(joinpath(dir, "$(cespgrn_prefix)_output*")))
    params = Dict(k => map(x -> map(x -> parse(Float64, x), split(x, "_")[[end-2, end-1, end]]), v) for (k, v) in files)
    all_nans = if static fill(NaN, Ng, Ng) else fill(NaN, Nc, Ng*Ng) end
    outputs = Dict( "locaTE_velo_dot"  => map(x -> try npzread(joinpath(x, "$(locate_what)_velo_dot.npy")) catch _ all_nans end, files["locaTE_velo_dot"]),
                   "locaTE_velo_cos"  => map(x -> try npzread(joinpath(x, "$(locate_what)_velo_cos.npy")) catch _ all_nans end, files["locaTE_velo_cos"]),
                   "locaTE_velo_corr" => map(x -> try npzread(joinpath(x, "$(locate_what)_velo_corr.npy")) catch _ all_nans end, files["locaTE_velo_corr"]),
                   "locaTE_dpt"       => map(x -> try npzread(joinpath(x, "$(locate_what)_dpt.npy")) catch _ all_nans end, files["locaTE_dpt"]),
                   "locaTE_statot"    => map(x -> try npzread(joinpath(x, "$(locate_what)_statot.npy")) catch _ all_nans end, files["locaTE_statot"]),
                   "locaTE_statot_ent"=> map(x -> try npzread(joinpath(x, "$(locate_what)_statot_ent.npy")) catch _ all_nans end, files["locaTE_statot_ent"]),
                   "locaTE_pba"       => map(x -> try npzread(joinpath(x, "$(locate_what)_pba.npy")) catch _ all_nans end, files["locaTE_pba"]),
                   "locaTE_undir"       => map(x -> try npzread(joinpath(x, "$(locate_what)_undir.npy")) catch _ all_nans end, files["locaTE_undir"]),
                   "cespgrn"         => map(x -> try preprocess_cespgrn(npzread(joinpath(x, "$(cespgrn_what)_cespgrn.npy"))) catch _  all_nans end, files["cespgrn"]))
    return Dict("files" => files, "params" => params, "outputs" => outputs)
end

function score_job_outputs(outputs, dataset; what = :aupr, kwargs...)
    out = Dict()
    for k in keys(outputs)
        ksymm = string(k, "_symm")
        if what == :aupr
            out[k] = map(x -> aupr(collect(eachcol(prec_rec_rate(dataset["R"]*dataset["J_mat"], x, Nq)))...; kwargs...), outputs[k])
            out[ksymm] = map(x -> aupr(collect(eachcol(prec_rec_rate(dataset["R"]*dataset["J_symm_mat"], locaTE.symm_row(x, Ng), Nq)))...; kwargs...), outputs[k])
        elseif what == :ep
            out[k] = map(x -> ep(collect(eachcol(prec_rec_rate(dataset["R"]*dataset["J_mat"], x, Nq)))...; kwargs...), outputs[k])
            out[ksymm] = map(x -> ep(collect(eachcol(prec_rec_rate(dataset["R"]*dataset["J_symm_mat"], locaTE.symm_row(x, Ng), Nq)))...; kwargs...), outputs[k])
        end
    end
    out
end

function score_map(score, data)
    function map(params, scores)
        p = hcat(params...)
        p_unique = unique.(eachrow(p))
        y = fill(Inf, length.(p_unique)...)
        for (i, x) in enumerate(eachcol(p))
            y[findall(s -> s == x[1], p_unique[1])[1],
              findall(s -> s == x[2], p_unique[2])[1],
              findall(s -> s == x[3], p_unique[3])[1]] = scores[i]
        end
        y
    end
    out = Dict()
    for k in keys(data["params"])
        ksymm = string(k, "_symm")
        out[k] = map(data["params"][k], score[k])
        out[ksymm] = map(data["params"][k], score[ksymm])
    end
    out
end

compute_baseline(x, thresh) = mean(x .> thresh)

function process_tenet(path, J; symm = false)
    G_tenet = CSV.read(path, DataFrame) 
    # deletecols!(G_tenet, 1)
    G_tenet = Array(G_tenet[:, 2:end])
    # G_tenet = locaTE.CLR(Array(G_tenet))
    # compute AUPR w.r.t ground truth
    G_tenet = symm ? locaTE.symm(G_tenet) : G_tenet
    p, r = collect(eachcol(prec_rec_rate(J, G_tenet, 512)))
    G_tenet, aupr(p, r)
end

function process_scode(path, J; symm = false)
    G_scode = abs.(Array(CSV.read(path, DataFrame; header = false)));
    G_scode = symm ? locaTE.symm(G_scode) : G_scode
    # compute AUPR w.r.t ground truth
    p, r = collect(eachcol(prec_rec_rate(J, G_scode, 512)))
    G_scode, aupr(p, r)
end

function process_locate(G, J; what = :aupr, symm = false, kwargs...)
    G_static = symm ? locaTE.symm(G) : G
    p, r = collect(eachcol(prec_rec_rate(J, G_static, 512)))
    if what == :aupr
        G_static, aupr(p, r)
    elseif what == :ep
        G_static, ep(p, r; kwargs...)
    end
end

function process_scribe(path, J; symm = false)
    G_scribe = Array(npzread(path));
    G_scribe[diagind(G_scribe)] .= 0
    G_scribe = locaTE.CLR(G_scribe);
    G_scribe = symm ? locaTE.symm(G_scribe) : G_scribe
    p, r = collect(eachcol(prec_rec_rate(J, G_scribe, 512)))
    G_scribe, aupr(p, r)
end

function process_pidc(path, J)
    G_pidc = Array(npzread(path));
    p, r = collect(eachcol(prec_rec_rate(J, G_pidc, 512)))
    G_pidc, aupr(p, r)
end

function process_grisli(path, J; symm = false)
    G_grisli = (-1)*abs.(Array(CSV.read(path, DataFrame; header = false)));
    G_grisli = symm ? locaTE.symm(G_grisli) : G_grisli
    # compute AUPR w.r.t ground truth
    p, r = collect(eachcol(prec_rec_rate(J, G_grisli, 512)))
    G_grisli, aupr(p, r)
end

function process_sincerities(path, J; symm = false)
    df_sincerities = CSV.read(path, DataFrame)
    G_sincerities = zeros(Ng, Ng)
    for (source, target, score, _) in eachrow(df_sincerities)
        i = parse(Int, split(source, "Gene ")[end])
        j = parse(Int, split(target, "Gene ")[end])
        G_sincerities[i, j] = score
    end
    G_sincerities = symm ? locaTE.symm(G_sincerities) : G_sincerities
    # compute AUPR w.r.t ground truth
    p, r = collect(eachcol(prec_rec_rate(J, G_sincerities, 512)))
    G_sincerities, aupr(p, r)
end

function process_genie3(path, J)
    G_genie3 = Array(CSV.read(path, DataFrame)[:, 2:end]);
    # compute AUPR w.r.t ground truth
    p, r = collect(eachcol(prec_rec_rate(J, G_genie3, 512)))
    G_genie3, aupr(p, r)
end
