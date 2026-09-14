using CSV
using DataFrames
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots
using Statistics: median

include(joinpath(@__DIR__, "AttractorSwitchingPOC.jl"))
using .AttractorSwitchingPOC

const CONDITION_SEED_STRIDE = 1_000_000

function wilson_interval(k::Integer, n::Integer; z::Float64 = 1.96)
    n > 0 || return (NaN, NaN)
    probability = k / n
    denominator = 1.0 + z^2 / n
    center = (probability + z^2 / (2n)) / denominator
    half_width =
        z * sqrt((probability * (1.0 - probability) + z^2 / (4n)) / n) /
        denominator
    return (
        max(0.0, center - half_width),
        min(1.0, center + half_width),
    )
end

function add_gate!(
    gates::DataFrame;
    gate::AbstractString,
    passed::Bool,
    metric::Real,
    target::Real,
    n::Integer,
    n_invalid::Integer,
    details::AbstractString,
)
    push!(
        gates,
        (
            gate = String(gate),
            passed = passed,
            metric = Float64(metric),
            target = Float64(target),
            n = Int(n),
            n_invalid = Int(n_invalid),
            details = String(details),
        ),
    )
    return gates
end

function parent_validation_ready(root::AbstractString)
    path = joinpath(
        root,
        "results",
        "post_meeting",
        "sde_validation",
        "validation_gates.csv",
    )
    isfile(path) || return (false, "missing validation_gates.csv")
    gates = CSV.read(path, DataFrame)
    required = Set([
        "zero_noise_matches_ode",
        "small_noise_ensemble_mean",
        "seed_reproducibility",
        "nonnegative_domain",
        "paper_regression",
        "time_resolution_sensitivity",
    ])
    observed = Set(String.(gates.gate))
    missing = setdiff(required, observed)
    isempty(missing) ||
        return (false, "missing parent gates: $(join(sort!(collect(missing)), ", "))")
    relevant = gates[in.(String.(gates.gate), Ref(required)), :]
    all(relevant.passed) || return (false, "a parent SDE gate failed")
    all(relevant.n_invalid .== 0) ||
        return (false, "parent SDE validation contains invalid paths")
    :source_fingerprint in propertynames(relevant) ||
        return (false, "parent validation lacks a source fingerprint")
    current_hash =
        AttractorSwitchingPOC.EffectorBSDE.source_fingerprint()
    all(String.(relevant.source_fingerprint) .== current_hash) ||
        return (false, "parent SDE source fingerprint is stale")
    :julia_version in propertynames(relevant) ||
        return (false, "parent validation lacks Julia provenance")
    parent_versions = unique(VersionNumber.(String.(relevant.julia_version)))
    same_minor_runtime = all(
        version ->
            version.major == VERSION.major && version.minor == VERSION.minor,
        parent_versions,
    )
    same_minor_runtime ||
        return (false, "parent validation used a different Julia major/minor runtime")
    parent_version_text = join(string.(parent_versions), ";")
    return (
        true,
        "all six parent SDE gates pass; source fingerprint matches; " *
        "parent Julia=$(parent_version_text), current Julia=$(VERSION)",
    )
end

function summarize_replicates(raw::DataFrame)
    summary = DataFrame()
    for group in groupby(raw, [:initial_attractor, :noise_scale])
        valid = group[group.valid .== true, :]
        n_requested = nrow(group)
        n_valid = nrow(valid)
        n_switched = count(valid.switched)
        lower, upper = wilson_interval(n_switched, n_valid)
        initial = first(group.initial_attractor)
        direction = initial == "low" ? "low_to_high" : "high_to_low"
        push!(
            summary,
            (
                direction = direction,
                initial_attractor = initial,
                target_attractor = initial == "low" ? "high" : "low",
                b6 = first(group.b6),
                noise_scale = first(group.noise_scale),
                noise_horizon = first(group.noise_horizon),
                relax_horizon = first(group.relax_horizon),
                n_requested = n_requested,
                n_valid = n_valid,
                n_invalid = n_requested - n_valid,
                n_switched = n_switched,
                p_switched = n_valid == 0 ? NaN : n_switched / n_valid,
                ci95_lower = lower,
                ci95_upper = upper,
                n_remained = n_valid - n_switched,
                median_nearest_distance = n_valid == 0 ? NaN :
                                          median(valid.nearest_distance),
            ),
        )
    end
    sort!(summary, [:direction, :noise_scale])
    return summary
end

function probability_figure(
    summary::DataFrame,
    config::SwitchingConfig,
    outpath::AbstractString,
)
    low = summary[summary.direction .== "low_to_high", :]
    high = summary[summary.direction .== "high_to_low", :]
    sort!(low, :noise_scale)
    sort!(high, :noise_scale)

    figure = plot(
        xlabel = "Noise amplitude epsilon (day^-1/2)",
        ylabel = "Switching probability",
        ylims = (0.0, 1.0),
        title = "Noise-induced attractor switching at b6 = $(config.b6)",
        legend = :topleft,
        size = (920, 610),
        left_margin = 7Plots.mm,
        bottom_margin = 7Plots.mm,
    )
    if nrow(low) > 0
        plot!(
            figure,
            low.noise_scale,
            low.p_switched;
            ribbon = (
                low.p_switched .- low.ci95_lower,
                low.ci95_upper .- low.p_switched,
            ),
            marker = :circle,
            markersize = 6,
            linewidth = 3,
            color = :darkorange,
            fillalpha = 0.18,
            label = "low tumor -> high tumor",
        )
    end
    if nrow(high) > 0
        plot!(
            figure,
            high.noise_scale,
            high.p_switched;
            ribbon = (
                high.p_switched .- high.ci95_lower,
                high.ci95_upper .- high.p_switched,
            ),
            marker = :diamond,
            markersize = 6,
            linewidth = 3,
            color = :steelblue,
            fillalpha = 0.18,
            label = "high tumor -> low tumor",
        )
    end
    annotate!(
        figure,
        maximum(summary.noise_scale) * 0.98,
        0.05,
        text("$(config.noise_horizon)-day noise pulse; deterministic relaxation afterward", 9, :right),
    )
    savefig(figure, outpath)
    return figure
end

function choose_representatives(raw::DataFrame)
    rows = raw[1:0, :]
    low = raw[raw.initial_attractor .== "low", :]
    high = raw[raw.initial_attractor .== "high", :]
    switching_levels = sort(unique(low[low.switched .== true, :noise_scale]))
    target_level = isempty(switching_levels) ? maximum(low.noise_scale) :
                   first(switching_levels)

    low_switch = low[
        (low.noise_scale .== target_level) .& low.switched .& low.valid,
        :,
    ]
    low_stay = low[
        (low.noise_scale .== target_level) .& .!low.switched .& low.valid,
        :,
    ]
    high_level = maximum(high.noise_scale)
    high_stay = high[
        (high.noise_scale .== high_level) .& .!high.switched .& high.valid,
        :,
    ]

    nrow(low_switch) > 0 && push!(rows, low_switch[1, :])
    nrow(low_stay) > 0 && push!(rows, low_stay[1, :])
    nrow(high_stay) > 0 && push!(rows, high_stay[1, :])
    return rows
end

function representative_figure(
    raw::DataFrame,
    attractors::AttractorSet,
    config::SwitchingConfig,
    outdir::AbstractString,
)
    selected = choose_representatives(raw)
    nrow(selected) > 0 || return DataFrame()
    panels = Any[]
    metadata = DataFrame()

    for row in eachrow(selected)
        initial = Symbol(row.initial_attractor)
        run = run_switching_once(
            initial,
            attractors;
            noise_scale = row.noise_scale,
            seed = row.seed,
            config = config,
        )
        title_text = run.switched ?
                     "$(run.initial_attractor) start -> $(run.final_attractor)" :
                     "$(run.initial_attractor) start -> remained $(run.final_attractor)"
        panel = plot(
            run.noise_solution.t,
            max.(run.noise_solution[1, :], 1.0e-8);
            yscale = :log10,
            color = :darkorange,
            linewidth = 2.5,
            label = "noise on",
            xlabel = "Time (days)",
            ylabel = "Tumor population",
            title = "$(title_text), epsilon=$(row.noise_scale)",
            legend = :right,
        )
        if !isnothing(run.relaxation_solution)
            plot!(
                panel,
                config.noise_horizon .+ run.relaxation_solution.t,
                max.(run.relaxation_solution[1, :], 1.0e-8);
                color = :steelblue,
                linewidth = 2.5,
                label = "noise off",
            )
        end
        hline!(
            panel,
            [attractors.low[1]];
            color = :forestgreen,
            linestyle = :dash,
            label = "low attractor",
        )
        hline!(
            panel,
            [attractors.saddle[1]];
            color = :black,
            linestyle = :dot,
            label = "saddle T (reference only)",
        )
        hline!(
            panel,
            [attractors.high[1]];
            color = :firebrick,
            linestyle = :dash,
            label = "high attractor",
        )
        vline!(
            panel,
            [config.noise_horizon];
            color = :gray40,
            linestyle = :dashdot,
            label = false,
        )
        push!(panels, panel)
        push!(
            metadata,
            (
                initial_attractor = run.initial_attractor,
                final_attractor = run.final_attractor,
                switched = run.switched,
                noise_scale = run.noise_scale,
                seed = run.seed,
            ),
        )
    end

    figure = plot(
        panels...;
        layout = (length(panels), 1),
        size = (980, 430 * length(panels)),
        left_margin = 7Plots.mm,
        bottom_margin = 6Plots.mm,
    )
    savefig(
        figure,
        joinpath(outdir, "representative_switching_trajectories.png"),
    )
    return metadata
end

function resolution_sensitivity(
    attractors::AttractorSet;
    noise_scale::Float64,
    n::Integer,
    base_seed::Integer,
    primary_config::SwitchingConfig,
)
    fine_config = SwitchingConfig(
        b6 = primary_config.b6,
        noise_horizon = primary_config.noise_horizon,
        relax_horizon = primary_config.relax_horizon,
        saveat = primary_config.saveat,
        relax_saveat = primary_config.relax_saveat,
        dt = 0.005,
        dtmax = 0.025,
        adaptive = true,
        reltol = primary_config.reltol,
        abstol = primary_config.abstol,
        maxiters = primary_config.maxiters,
        classification_tolerance = primary_config.classification_tolerance,
    )
    rows = DataFrame()
    for replicate in 1:n
        seed = base_seed + replicate
        primary = run_switching_once(
            :low,
            attractors;
            noise_scale = noise_scale,
            seed = seed,
            config = primary_config,
        )
        fine = run_switching_once(
            :low,
            attractors;
            noise_scale = noise_scale,
            seed = seed,
            config = fine_config,
        )
        push!(
            rows,
            (
                replicate = replicate,
                seed = seed,
                noise_scale = noise_scale,
                primary_valid = primary.valid,
                fine_valid = fine.valid,
                primary_switched = primary.switched,
                fine_switched = fine.switched,
                classifications_agree =
                    primary.valid && fine.valid &&
                    primary.final_attractor == fine.final_attractor,
            ),
        )
    end
    valid = rows[rows.primary_valid .& rows.fine_valid, :]
    primary_probability = isempty(valid.primary_switched) ? NaN :
                          sum(valid.primary_switched) / nrow(valid)
    fine_probability = isempty(valid.fine_switched) ? NaN :
                       sum(valid.fine_switched) / nrow(valid)
    difference = abs(primary_probability - fine_probability)
    agreement = nrow(valid) == 0 ? NaN :
                sum(valid.classifications_agree) / nrow(valid)
    passed =
        nrow(valid) == n && isfinite(difference) && difference <= 0.10 &&
        agreement >= 0.90
    return (
        rows = rows,
        passed = passed,
        probability_difference = difference,
        agreement = agreement,
        n_invalid = n - nrow(valid),
    )
end

function write_text_summary(
    path::AbstractString,
    summary::DataFrame,
    gates::DataFrame,
    config::SwitchingConfig,
    source_hash::AbstractString,
    git_sha::AbstractString,
)
    open(path, "w") do io
        println(io, "Attractor-switching proof-of-concept")
        println(io, "=====================================")
        println(io)
        println(io, "b6: ", config.b6)
        println(io, "Noise pulse: ", config.noise_horizon, " days")
        println(io, "Deterministic relaxation: ", config.relax_horizon, " days")
        println(io, "Noise model: common multiplicative Ito noise")
        println(io, "Solver: SOSRI")
        println(io, "Git SHA: ", git_sha)
        println(io, "POC source fingerprint: ", source_hash)
        println(io, "Julia: ", VERSION)
        println(io)
        println(io, "Validation gates: ", count(gates.passed), "/", nrow(gates), " passed")
        println(io)
        println(io, "Switching estimates:")
        for row in eachrow(summary)
            println(
                io,
                "  ",
                row.direction,
                ", epsilon=",
                row.noise_scale,
                ": ",
                row.n_switched,
                "/",
                row.n_valid,
                " = ",
                round(row.p_switched; digits = 3),
                " (95% CI ",
                round(row.ci95_lower; digits = 3),
                "-",
                round(row.ci95_upper; digits = 3),
                ")",
            )
        end
        println(io)
        println(io, "Interpretation boundary:")
        println(io, "  This is a numerical proof-of-concept, not a calibrated biological")
        println(io, "  estimate of in-vivo switching frequency or noise amplitude.")
    end
end

function run_poc(;
    nsims::Integer = 100,
    resolution_reps::Integer = 50,
    base_seed::Integer = 20260805,
    noise_scales::Vector{Float64} = [0.0, 0.10, 0.15, 0.20, 0.25, 0.30],
    config::SwitchingConfig = SwitchingConfig(),
    outdir::AbstractString = joinpath(
        @__DIR__,
        "..",
        "..",
        "..",
        "results",
        "post_meeting",
        "attractor_switching_poc",
    ),
)
    nsims >= 50 || throw(ArgumentError("nsims must be at least 50"))
    resolution_reps >= 30 ||
        throw(ArgumentError("resolution_reps must be at least 30"))
    all(noise_scales .>= 0.0) ||
        throw(ArgumentError("noise scales must be nonnegative"))
    0.0 in noise_scales ||
        throw(ArgumentError("noise grid must include epsilon=0"))
    mkpath(outdir)
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    continuation_path = joinpath(
        root,
        "results",
        "post_meeting",
        "formal_b6_continuation",
        "formal_continuation_all_points.csv",
    )
    attractors = refine_attractor_set(continuation_path; b6 = config.b6)
    source_hash = poc_source_fingerprint()
    git_sha = repository_git_sha()
    CSV.write(
        joinpath(outdir, "equilibria.csv"),
        equilibrium_table(attractors, config.b6),
    )

    gates = DataFrame()
    parent_ready, parent_details = parent_validation_ready(root)
    add_gate!(
        gates;
        gate = "parent_sde_validation_current",
        passed = parent_ready,
        metric = parent_ready ? 1.0 : 0.0,
        target = 1.0,
        n = 6,
        n_invalid = 0,
        details = parent_details,
    )

    equilibrium_ok =
        maximum([
            attractors.low_residual,
            attractors.saddle_residual,
            attractors.high_residual,
        ]) < 1.0e-9 &&
        attractors.low_n_unstable == 0 &&
        attractors.saddle_n_unstable == 1 &&
        attractors.high_n_unstable == 0
    add_gate!(
        gates;
        gate = "three_equilibria_and_stability",
        passed = equilibrium_ok,
        metric = maximum([
            attractors.low_residual,
            attractors.saddle_residual,
            attractors.high_residual,
        ]),
        target = 1.0e-9,
        n = 3,
        n_invalid = 0,
        details = "Low and high equilibria must be stable; the middle equilibrium must have exactly one unstable mode.",
    )

    zero_runs = [
        run_switching_once(
            initial,
            attractors;
            noise_scale = 0.0,
            seed = base_seed + index,
            config = config,
        ) for (index, initial) in enumerate((:low, :high))
    ]
    zero_ok = all(run -> run.valid && !run.switched, zero_runs)
    zero_error = maximum(
        run.classification.nearest_distance for run in zero_runs if
        !isnothing(run.classification)
    )
    add_gate!(
        gates;
        gate = "zero_noise_preserves_attractors",
        passed = zero_ok && zero_error < config.classification_tolerance,
        metric = zero_error,
        target = config.classification_tolerance,
        n = 2,
        n_invalid = count(run -> !run.valid, zero_runs),
        details = "At epsilon=0, both equilibria must remain in their original deterministic basins.",
    )

    reproducibility_1 = solve_noise_pulse(
        attractors.low;
        noise_scale = 0.20,
        seed = base_seed + 10_000,
        config = config,
    )
    reproducibility_2 = solve_noise_pulse(
        attractors.low;
        noise_scale = 0.20,
        seed = base_seed + 10_000,
        config = config,
    )
    reproducible =
        reproducibility_1.t == reproducibility_2.t &&
        reproducibility_1.u == reproducibility_2.u
    add_gate!(
        gates;
        gate = "seed_reproducibility",
        passed = reproducible,
        metric = reproducible ? 0.0 : 1.0,
        target = 0.0,
        n = 2,
        n_invalid = 0,
        details = "The same explicit seed must reproduce the saved stochastic path exactly.",
    )

    preflight_path = joinpath(outdir, "validation_gates.csv")
    preflight = copy(gates)
    preflight.git_sha = fill(git_sha, nrow(preflight))
    preflight.source_fingerprint = fill(source_hash, nrow(preflight))
    preflight.julia_version = fill(string(VERSION), nrow(preflight))
    CSV.write(preflight_path, preflight)
    all(gates.passed) ||
        error("Preflight validation failed; switching ensemble is blocked")

    raw = DataFrame()
    run_id = 0
    conditions = collect(Iterators.product((:low, :high), noise_scales))
    for (condition_index, (initial, noise_scale)) in enumerate(conditions)
        println(
            "Running start=",
            initial,
            ", epsilon=",
            noise_scale,
            ", nsims=",
            nsims,
        )
        for replicate in 1:nsims
            run_id += 1
            seed =
                base_seed + condition_index * CONDITION_SEED_STRIDE + replicate
            run = run_switching_once(
                initial,
                attractors;
                noise_scale = noise_scale,
                seed = seed,
                config = config,
            )
            push!(
                raw,
                merge(
                    (
                        run_id = run_id,
                        replicate = replicate,
                        git_sha = git_sha,
                        source_fingerprint = source_hash,
                        julia_version = string(VERSION),
                        solver = "SOSRI",
                        noise_mode = "paper_common",
                        gompertz_mode = "continuation",
                        dt_initial = config.dt,
                        dtmax = config.dtmax,
                        adaptive = config.adaptive,
                    ),
                    run_to_metrics(run, config),
                ),
            )
        end
    end
    raw_path = joinpath(outdir, "switching_replicates.csv")
    CSV.write(raw_path, raw)
    summary = summarize_replicates(raw)
    CSV.write(joinpath(outdir, "switching_summary.csv"), summary)

    invalid_count = count(.!raw.valid)
    add_gate!(
        gates;
        gate = "nonnegative_full_horizon_paths",
        passed = invalid_count == 0 && minimum(raw.min_state_saved) >= 0.0,
        metric = minimum(raw.min_state_saved),
        target = 0.0,
        n = nrow(raw),
        n_invalid = invalid_count,
        details = "Every SDE path must remain finite and nonnegative and reach the full noise horizon.",
    )
    ambiguous_count = count(raw.final_attractor .== "ambiguous")
    add_gate!(
        gates;
        gate = "deterministic_basin_classification",
        passed = ambiguous_count == 0,
        metric = maximum(raw.nearest_distance),
        target = config.classification_tolerance,
        n = nrow(raw),
        n_invalid = ambiguous_count,
        details = "After noise is removed, every valid endpoint must relax within tolerance of one deterministic attractor.",
    )

    sensitivity = resolution_sensitivity(
        attractors;
        noise_scale = 0.20,
        n = resolution_reps,
        base_seed = base_seed + 50_000_000,
        primary_config = config,
    )
    CSV.write(
        joinpath(outdir, "time_resolution_sensitivity.csv"),
        sensitivity.rows,
    )
    add_gate!(
        gates;
        gate = "switching_time_resolution_sensitivity",
        passed = sensitivity.passed,
        metric = sensitivity.probability_difference,
        target = 0.10,
        n = resolution_reps,
        n_invalid = sensitivity.n_invalid,
        details = "At epsilon=0.20, primary and finer adaptive settings must differ by <=0.10 with >=90% per-seed fate agreement (observed $(round(sensitivity.agreement; digits=3))).",
    )

    gates.git_sha = fill(git_sha, nrow(gates))
    gates.source_fingerprint = fill(source_hash, nrow(gates))
    gates.julia_version = fill(string(VERSION), nrow(gates))
    CSV.write(preflight_path, gates)

    probability_figure(
        summary,
        config,
        joinpath(outdir, "attractor_switching_probability.png"),
    )
    representative_metadata =
        representative_figure(raw, attractors, config, outdir)
    CSV.write(
        joinpath(outdir, "representative_trajectories.csv"),
        representative_metadata,
    )

    configuration = DataFrame(
        setting = [
            "b6",
            "noise_horizon_days",
            "relax_horizon_days",
            "noise_scales",
            "replicates_per_condition",
            "resolution_replicates",
            "base_seed",
            "noise_mode",
            "solver",
            "dt_initial",
            "dtmax",
            "classification_tolerance",
        ],
        value = string.([
            config.b6,
            config.noise_horizon,
            config.relax_horizon,
            join(noise_scales, ";"),
            nsims,
            resolution_reps,
            base_seed,
            "paper_common",
            "SOSRI",
            config.dt,
            config.dtmax,
            config.classification_tolerance,
        ]),
    )
    CSV.write(joinpath(outdir, "experiment_configuration.csv"), configuration)
    write_text_summary(
        joinpath(outdir, "poc_summary.txt"),
        summary,
        gates,
        config,
        source_hash,
        git_sha,
    )

    all(gates.passed) ||
        error("At least one final validation gate failed; inspect validation_gates.csv")
    println("All ", nrow(gates), " attractor-switching validation gates passed.")
    println("Outputs written to ", outdir)
    return (
        attractors = attractors,
        raw = raw,
        summary = summary,
        gates = gates,
    )
end

function parse_noise_scales(value::AbstractString)
    return parse.(Float64, strip.(split(value, ",")))
end

if abspath(PROGRAM_FILE) == @__FILE__
    nsims = parse(Int, get(ENV, "POC_NSIMS", "100"))
    resolution_reps =
        parse(Int, get(ENV, "POC_RESOLUTION_REPS", "50"))
    noise_scales = parse_noise_scales(
        get(ENV, "POC_NOISE_SCALES", "0.0,0.10,0.15,0.20,0.25,0.30"),
    )
    run_poc(
        nsims = nsims,
        resolution_reps = resolution_reps,
        noise_scales = noise_scales,
    )
end
