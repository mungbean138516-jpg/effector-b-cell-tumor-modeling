using CSV
using DataFrames
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots
using Statistics

include(joinpath(@__DIR__, "EffectorBSDE.jl"))
using .EffectorBSDE

function wilson_interval(k::Integer, n::Integer; z = 1.96)
    n > 0 || return (NaN, NaN)
    p = k / n
    denominator = 1 + z^2 / n
    center = (p + z^2 / (2n)) / denominator
    half =
        z * sqrt((p * (1 - p) + z^2 / (4n)) / n) / denominator
    return (max(0.0, center - half), min(1.0, center + half))
end

function assert_validation_ready(
    outdir::AbstractString,
    source_hash::AbstractString,
)
    validation_path = joinpath(outdir, "validation_gates.csv")
    isfile(validation_path) ||
        error(
            "Missing validation_gates.csv. Run validate_sde_pipeline.jl " *
            "before the pilot.",
        )
    gates = CSV.read(validation_path, DataFrame)
    required = Set([
        "zero_noise_matches_ode",
        "small_noise_ensemble_mean",
        "seed_reproducibility",
        "nonnegative_domain",
        "paper_regression",
        "time_resolution_sensitivity",
    ])
    observed = Set(gates.gate)
    missing = setdiff(required, observed)
    if !isempty(missing)
        missing_text = join(sort!(collect(missing)), ", ")
        error("Validation artifact is missing gates: $(missing_text).")
    end
    relevant = gates[in.(gates.gate, Ref(required)), :]
    all(relevant.passed) ||
        error("One or more validation gates failed; pilot is blocked.")
    all(relevant.n_invalid .== 0) ||
        error("Validation contains invalid trajectories; pilot is blocked.")
    :source_fingerprint in propertynames(relevant) ||
        error("Validation artifact lacks a source fingerprint.")
    all(string.(relevant.source_fingerprint) .== source_hash) ||
        error(
            "Validation artifact does not match the current model, solver, " *
            "or dependency files. Rerun validation before the pilot.",
        )
    :julia_version in propertynames(relevant) ||
        error("Validation artifact lacks Julia runtime provenance.")
    all(string.(relevant.julia_version) .== string(VERSION)) ||
        error(
            "Validation used a different Julia runtime. Rerun validation " *
            "under Julia $(VERSION) before the pilot.",
        )
    return validation_path
end

function summarize_replicates(raw::DataFrame)
    rows = DataFrame()
    for group in groupby(raw, [:b6_multiplier, :b6])
        valid = group[group.valid .== true, :]
        n_requested = nrow(group)
        n_valid = nrow(valid)
        n_paper = count(valid.paper_established)
        n_strict = count(valid.strict_established)
        n_intermediate = count(valid.strict_state .== "persistent_intermediate")
        paper_lo, paper_hi = wilson_interval(n_paper, n_valid)
        strict_lo, strict_hi = wilson_interval(n_strict, n_valid)
        established = valid[valid.paper_established .== true, :]
        removed_times = collect(skipmissing(valid.removal_time))
        removed_times = filter(isfinite, removed_times)

        push!(
            rows,
            (
                b6_multiplier = first(group.b6_multiplier),
                b6 = first(group.b6),
                n_requested = n_requested,
                n_valid = n_valid,
                n_invalid = n_requested - n_valid,
                n_paper_established = n_paper,
                p_paper_established =
                    n_valid == 0 ? NaN : n_paper / n_valid,
                p_paper_lo = paper_lo,
                p_paper_hi = paper_hi,
                n_strict_established = n_strict,
                n_persistent_intermediate = n_intermediate,
                p_strict_established =
                    n_valid == 0 ? NaN : n_strict / n_valid,
                p_strict_lo = strict_lo,
                p_strict_hi = strict_hi,
                median_T_final_all_valid =
                    n_valid == 0 ? NaN : median(valid.T_final),
                median_T_final_paper_established =
                    isempty(established.T_final) ? NaN :
                    median(established.T_final),
                median_T_max_all_valid =
                    n_valid == 0 ? NaN : median(valid.T_max_saved),
                median_removal_time =
                    isempty(removed_times) ? NaN : median(removed_times),
                total_solver_rejections =
                    n_valid == 0 ? 0 : sum(valid.solver_rejections),
            ),
        )
    end
    sort!(rows, :b6_multiplier)
    return rows
end

function build_summary_figure(summary::DataFrame, outpath::AbstractString)
    x = summary.b6_multiplier
    paper_lower = summary.p_paper_established .- summary.p_paper_lo
    paper_upper = summary.p_paper_hi .- summary.p_paper_established
    strict_lower = summary.p_strict_established .- summary.p_strict_lo
    strict_upper = summary.p_strict_hi .- summary.p_strict_established

    p1 = plot(
        x,
        summary.p_paper_established;
        ribbon = (paper_lower, paper_upper),
        marker = :circle,
        linewidth = 2.5,
        label = "CIR paper criterion",
        xlabel = "b6 / b6*",
        ylabel = "Probability",
        ylims = (0, 1),
        title = "One-year stochastic establishment",
        legend = :best,
    )
    plot!(
        p1,
        x,
        summary.p_strict_established;
        ribbon = (strict_lower, strict_upper),
        marker = :diamond,
        linewidth = 2.5,
        label = "Strict sensitivity criterion",
    )
    vline!(p1, [1.0]; linestyle = :dash, color = :black, label = "b6*")

    p2 = plot(
        x,
        max.(summary.median_T_final_paper_established, 1.0e-6);
        marker = :circle,
        linewidth = 2.5,
        yscale = :log10,
        xlabel = "b6 / b6*",
        ylabel = "Median final tumor",
        title = "Final burden among paper-established trajectories",
        label = false,
    )
    vline!(p2, [1.0]; linestyle = :dash, color = :black, label = false)

    fig = plot(p1, p2; layout = (2, 1), size = (900, 850))
    savefig(fig, outpath)
    return fig
end

function run_near_fold_pilot(;
    nsims::Integer = 100,
    base_seed::Integer = 20260728,
    b5::Float64 = 1.0e-4,
    noise_scale::Float64 = 1.0,
    dtmax::Float64 = 0.05,
    outdir::AbstractString = joinpath(
        @__DIR__,
        "..",
        "..",
        "..",
        "results",
        "post_meeting",
        "sde_validation",
    ),
    git_sha::AbstractString = repository_git_sha(),
    source_hash::AbstractString = source_fingerprint(),
    require_validation::Bool = true,
)
    nsims >= 100 || throw(ArgumentError("nsims must be at least 100"))
    b5 == 1.0e-4 ||
        throw(ArgumentError("validated pilot requires b5=1e-4"))
    noise_scale == 1.0 ||
        throw(ArgumentError("validated pilot requires noise_scale=1.0"))
    dtmax == 0.05 ||
        throw(ArgumentError("validated pilot requires dtmax=0.05"))
    mkpath(outdir)
    require_validation && assert_validation_ready(outdir, source_hash)

    multipliers = [0.8, 0.9, 1.0, 1.1, 1.2]
    b6_values = DEFAULT_B6_STAR .* multipliers
    config = SimulationConfig(
        horizon = 365.0,
        saveat = 1.0,
        dt = min(0.01, dtmax),
        dtmax = dtmax,
        noise_scale = noise_scale,
        noise_mode = :paper_common,
        gompertz_mode = :continuation,
        paper_threshold = 1.0,
        strict_threshold = 1.0e3,
        strict_window_days = 30.0,
        strict_required_fraction = 0.90,
    )

    raw = DataFrame()
    for (condition_index, (multiplier, b6)) in
        enumerate(zip(multipliers, b6_values))
        println(
            "Running b6/b6* = ",
            multiplier,
            " (b6 = ",
            b6,
            "), nsims = ",
            nsims,
        )
        for replicate in 1:nsims
            seed = base_seed + (condition_index - 1) * nsims + replicate
            run = solve_sde_once(
                b5 = b5,
                b6 = b6,
                seed = seed,
                config = config,
            )
            metrics = trajectory_metrics(run)
            push!(
                raw,
                merge(
                    (
                        run_id = (condition_index - 1) * nsims + replicate,
                        replicate = replicate,
                        git_sha = git_sha,
                        source_fingerprint = source_hash,
                        julia_version = string(VERSION),
                        b5 = b5,
                        b6 = b6,
                        b6_multiplier = multiplier,
                        b6_star = DEFAULT_B6_STAR,
                        noise_scale = noise_scale,
                        noise_mode = string(config.noise_mode),
                        gompertz_mode = string(config.gompertz_mode),
                        solver = "SOSRI",
                        dt_initial = config.dt,
                        dtmax = config.dtmax,
                        adaptive = config.adaptive,
                        maxiters = config.maxiters,
                        horizon = config.horizon,
                        paper_threshold = config.paper_threshold,
                        strict_threshold = config.strict_threshold,
                        strict_window_days = config.strict_window_days,
                        strict_required_fraction =
                            config.strict_required_fraction,
                    ),
                    metrics,
                ),
            )
        end
    end

    summary = summarize_replicates(raw)
    raw_path = joinpath(outdir, "near_fold_pilot_replicates.csv")
    summary_path = joinpath(outdir, "near_fold_pilot_summary.csv")
    figure_path = joinpath(outdir, "near_fold_establishment_probability.png")
    CSV.write(raw_path, raw)
    CSV.write(summary_path, summary)
    build_summary_figure(summary, figure_path)

    n_invalid = count(.!raw.valid)
    println("Saved: ", raw_path)
    println("Saved: ", summary_path)
    println("Saved: ", figure_path)
    println(summary)
    n_invalid == 0 ||
        error(
            "Pilot produced $(n_invalid) invalid trajectories; " *
            "do not interpret biological probabilities.",
        )
    return raw, summary
end

function main()
    nsims = parse(Int, get(ENV, "SDE_NSIMS", "100"))
    base_seed = parse(Int, get(ENV, "SDE_BASE_SEED", "20260728"))
    noise_scale = parse(Float64, get(ENV, "SDE_NOISE_SCALE", "1.0"))
    dtmax = parse(Float64, get(ENV, "SDE_DTMAX", "0.05"))
    return run_near_fold_pilot(
        nsims = nsims,
        base_seed = base_seed,
        noise_scale = noise_scale,
        dtmax = dtmax,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
