using CSV
using DataFrames
using Statistics

include(joinpath(@__DIR__, "EffectorBSDE.jl"))
using .EffectorBSDE

function wilson_interval(k::Integer, n::Integer; z = 2.575829)
    n > 0 || return (NaN, NaN)
    p = k / n
    denominator = 1 + z^2 / n
    center = (p + z^2 / (2n)) / denominator
    half =
        z * sqrt((p * (1 - p) + z^2 / (4n)) / n) / denominator
    return (max(0.0, center - half), min(1.0, center + half))
end

function normalized_max_error(sol_a, sol_b)
    size(sol_a) == size(sol_b) ||
        throw(ArgumentError("solutions must have matching saved grids"))
    errors = Float64[]
    for i in axes(sol_a, 1)
        scale = max(maximum(abs.(sol_b[i, :])), 1.0)
        push!(errors, maximum(abs.(sol_a[i, :] .- sol_b[i, :])) / scale)
    end
    return maximum(errors)
end

function paper_regression(;
    nsims::Integer,
    base_seed::Integer,
    dtmax::Float64,
    git_sha::AbstractString,
)
    config = SimulationConfig(
        horizon = 365.0,
        saveat = 1.0,
        dt = min(0.01, dtmax),
        dtmax = dtmax,
        noise_scale = 1.0,
        noise_mode = :paper_common,
        gompertz_mode = :continuation,
    )
    rows = DataFrame()
    for replicate in 1:nsims
        run = solve_sde_once(
            a10 = 0.0,
            a11 = 0.0,
            b6 = 0.0,
            seed = base_seed + replicate,
            config = config,
        )
        push!(
            rows,
            merge(
                (
                    replicate = replicate,
                    git_sha = git_sha,
                    julia_version = string(VERSION),
                    solver = "SOSRI",
                    noise_scale = config.noise_scale,
                    noise_mode = string(config.noise_mode),
                    gompertz_mode = string(config.gompertz_mode),
                    dt_initial = config.dt,
                    dtmax = config.dtmax,
                    adaptive = config.adaptive,
                    maxiters = config.maxiters,
                    horizon = config.horizon,
                ),
                trajectory_metrics(run),
            ),
        )
    end
    valid = rows[rows.valid .== true, :]
    successes = count(valid.paper_established)
    lo, hi = wilson_interval(successes, nrow(valid))
    reference = 0.321711
    passed =
        nrow(valid) == nsims && lo <= reference <= hi &&
        count(valid.negativity_violation) == 0
    return (
        rows = rows,
        result = (
            gate = "paper_regression",
            passed = passed,
            metric = nrow(valid) == 0 ? NaN : successes / nrow(valid),
            lower = lo,
            upper = hi,
            target = reference,
            n = nsims,
            n_invalid = nsims - nrow(valid),
            details =
                "Original 4-state effects recovered by setting a10=a11=b6=0; " *
                "99% Wilson interval should contain 0.321711.",
        ),
    )
end

function resolution_check(;
    nsims::Integer,
    base_seed::Integer,
    dt_coarse::Float64 = 0.01,
    dt_fine::Float64 = 0.005,
)
    probabilities = Float64[]
    invalid_counts = Int[]
    configs = (
        SimulationConfig(
            horizon = 365.0,
            saveat = 1.0,
            dt = dt_coarse,
            dtmax = dt_coarse,
            adaptive = false,
            noise_scale = 1.0,
            noise_mode = :paper_common,
            gompertz_mode = :continuation,
        ),
        SimulationConfig(
            horizon = 365.0,
            saveat = 1.0,
            dt = dt_fine,
            dtmax = dt_fine,
            adaptive = false,
            noise_scale = 1.0,
            noise_mode = :paper_common,
            gompertz_mode = :continuation,
        ),
        SimulationConfig(
            horizon = 365.0,
            saveat = 1.0,
            dt = dt_fine,
            dtmax = 0.05,
            adaptive = true,
            noise_scale = 1.0,
            noise_mode = :paper_common,
            gompertz_mode = :continuation,
        ),
    )
    for config in configs
        metrics = [
            trajectory_metrics(
                solve_sde_once(
                    b6 = 1.1 * DEFAULT_B6_STAR,
                    seed = base_seed + replicate,
                    config = config,
                ),
            ) for replicate in 1:nsims
        ]
        valid = filter(m -> m.valid, metrics)
        push!(
            probabilities,
            isempty(valid) ? NaN :
            count(m -> m.paper_established, valid) / length(valid),
        )
        push!(invalid_counts, nsims - length(valid))
    end
    fixed_difference = abs(probabilities[1] - probabilities[2])
    adaptive_difference = abs(probabilities[2] - probabilities[3])
    difference = max(fixed_difference, adaptive_difference)
    passed =
        all(==(0), invalid_counts) && isfinite(difference) &&
        difference <= 0.10
    return (
        gate = "time_resolution_sensitivity",
        passed = passed,
        metric = difference,
        lower = probabilities[1],
        upper = probabilities[2],
        target = 0.10,
        n = 3nsims,
        n_invalid = sum(invalid_counts),
        details =
            "Maximum probability difference across fixed dt=0.01, fixed " *
            "dt=0.005, and the primary adaptive solver. lower/upper fields " *
            "contain the two fixed-step estimates.",
    )
end

function run_validation_suite(;
    paper_reps::Integer = 500,
    resolution_reps::Integer = 100,
    base_seed::Integer = 20260728,
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
)
    paper_reps >= 500 ||
        throw(ArgumentError("paper_reps must be at least 500"))
    resolution_reps >= 100 ||
        throw(ArgumentError("resolution_reps must be at least 100"))
    mkpath(outdir)
    gates = DataFrame()

    zero_config = SimulationConfig(
        horizon = 60.0,
        saveat = 0.5,
        dt = 0.01,
        dtmax = 0.05,
        noise_scale = 0.0,
        noise_mode = :paper_common,
        gompertz_mode = :continuation,
    )
    zero_run_1 = solve_sde_once(
        b6 = 0.8 * DEFAULT_B6_STAR,
        seed = base_seed + 1,
        config = zero_config,
    )
    zero_run_2 = solve_sde_once(
        b6 = 0.8 * DEFAULT_B6_STAR,
        seed = base_seed + 2,
        config = zero_config,
    )
    ode = solve_ode_reference(
        b6 = 0.8 * DEFAULT_B6_STAR,
        config = zero_config,
    )
    zero_error = normalized_max_error(zero_run_1.solution, ode)
    seed_invariant =
        zero_run_1.solution.t == zero_run_2.solution.t &&
        zero_run_1.solution.u == zero_run_2.solution.u
    zero_metrics = trajectory_metrics(zero_run_1)
    push!(
        gates,
        (
            gate = "zero_noise_matches_ode",
            passed =
                zero_error < 1.0e-3 && seed_invariant &&
                zero_metrics.valid,
            metric = zero_error,
            lower = NaN,
            upper = NaN,
            target = 1.0e-3,
            n = 2,
            n_invalid = zero_metrics.valid ? 0 : 1,
            details =
                "Normalized max state error; zero-noise paths must also be " *
                "seed-invariant.",
        ),
    )

    small_noise_config = SimulationConfig(
        horizon = 30.0,
        saveat = 1.0,
        dt = 0.01,
        dtmax = 0.05,
        noise_scale = 0.025,
        noise_mode = :paper_common,
        gompertz_mode = :continuation,
    )
    small_noise_runs = [
        solve_sde_once(
            b6 = 0.8 * DEFAULT_B6_STAR,
            seed = base_seed + 500 + i,
            config = small_noise_config,
        ) for i in 1:50
    ]
    small_noise_metrics = trajectory_metrics.(small_noise_runs)
    small_noise_invalid = count(m -> !m.valid, small_noise_metrics)
    small_noise_mean = vec(
        mean(
            reduce(
                hcat,
                [run.solution.u[end] for run in small_noise_runs],
            );
            dims = 2,
        ),
    )
    small_noise_ode = solve_ode_reference(
        b6 = 0.8 * DEFAULT_B6_STAR,
        config = small_noise_config,
    )
    small_noise_reference = small_noise_ode.u[end]
    small_noise_error = maximum(
        abs.(small_noise_mean .- small_noise_reference) ./
        max.(abs.(small_noise_reference), 1.0),
    )
    push!(
        gates,
        (
            gate = "small_noise_ensemble_mean",
            passed =
                small_noise_invalid == 0 && small_noise_error < 0.10,
            metric = small_noise_error,
            lower = NaN,
            upper = NaN,
            target = 0.10,
            n = length(small_noise_runs),
            n_invalid = small_noise_invalid,
            details =
                "At noise_scale=0.025 and b6=0.8*b6*, normalized terminal " *
                "ensemble-mean error versus the ODE must remain below 0.10.",
        ),
    )

    reproducibility_config = SimulationConfig(
        horizon = 30.0,
        saveat = 0.5,
        dt = 0.01,
        dtmax = 0.05,
        noise_scale = 0.1,
        noise_mode = :paper_common,
        gompertz_mode = :continuation,
    )
    reproducible_1 = solve_sde_once(
        seed = base_seed + 10,
        config = reproducibility_config,
    )
    reproducible_2 = solve_sde_once(
        seed = base_seed + 10,
        config = reproducibility_config,
    )
    exact_repeat =
        reproducible_1.solution.t == reproducible_2.solution.t &&
        reproducible_1.solution.u == reproducible_2.solution.u
    reproducible_metrics = trajectory_metrics(reproducible_1)
    push!(
        gates,
        (
            gate = "seed_reproducibility",
            passed = exact_repeat && reproducible_metrics.valid,
            metric = exact_repeat ? 0.0 : 1.0,
            lower = NaN,
            upper = NaN,
            target = 0.0,
            n = 2,
            n_invalid = reproducible_metrics.valid ? 0 : 1,
            details = "Identical seed must reproduce saved times and states.",
        ),
    )

    nonnegative_metrics = [
        trajectory_metrics(
            solve_sde_once(
                seed = base_seed + 100 + i,
                config = SimulationConfig(
                    horizon = 30.0,
                    saveat = 0.5,
                    dt = 0.01,
                    dtmax = 0.05,
                    noise_scale = 1.0,
                    noise_mode = :paper_common,
                    gompertz_mode = :continuation,
                ),
            ),
        ) for i in 1:20
    ]
    nonnegative_failures = count(
        m -> !m.valid || m.negativity_violation,
        nonnegative_metrics,
    )
    push!(
        gates,
        (
            gate = "nonnegative_domain",
            passed = nonnegative_failures == 0,
            metric = nonnegative_failures,
            lower = NaN,
            upper = NaN,
            target = 0.0,
            n = length(nonnegative_metrics),
            n_invalid = nonnegative_failures,
            details =
                "No accepted saved state may be negative or nonfinite.",
        ),
    )

    regression = paper_regression(
        nsims = paper_reps,
        base_seed = base_seed + 1_000,
        dtmax = 0.05,
        git_sha = git_sha,
    )
    append!(gates, DataFrame([regression.result]))
    resolution = resolution_check(
        nsims = resolution_reps,
        base_seed = base_seed + 10_000,
    )
    append!(gates, DataFrame([resolution]))
    gates.git_sha = fill(git_sha, nrow(gates))
    gates.julia_version = fill(string(VERSION), nrow(gates))

    gate_path = joinpath(outdir, "validation_gates.csv")
    regression_path = joinpath(outdir, "paper_regression_replicates.csv")
    CSV.write(gate_path, gates)
    CSV.write(regression_path, regression.rows)
    println(gates)
    println("Saved: ", gate_path)
    println("Saved: ", regression_path)

    all(gates.passed) ||
        error("One or more SDE validation gates failed; do not run the pilot.")
    return gates, regression.rows
end

function main()
    paper_reps = parse(Int, get(ENV, "SDE_PAPER_REPS", "500"))
    resolution_reps =
        parse(Int, get(ENV, "SDE_RESOLUTION_REPS", "100"))
    return run_validation_suite(
        paper_reps = paper_reps,
        resolution_reps = resolution_reps,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
