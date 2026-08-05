module EffectorBSDE

using OrdinaryDiffEq: Tsit5
using SHA: sha256
using SciMLBase:
    CallbackSet,
    DiscreteCallback,
    ODEProblem,
    SDEProblem,
    solve,
    successful_retcode,
    terminate!
using StochasticDiffEq: SOSRI, WienerProcess

export DEFAULT_B6_STAR,
       SimulationConfig,
       base_parameters,
       classify_tumor_path,
       inoculation_start,
       model_B!,
       repository_git_sha,
       source_fingerprint,
       solve_ode_reference,
       solve_sde_once,
       trajectory_metrics

"""
Fold location from the repository's zero-eigenvalue classification.

This is a numerical estimate, not an exact biological constant.
"""
const DEFAULT_B6_STAR = 6.676254353545294e-7

function repository_git_sha()
    haskey(ENV, "GIT_SHA") && return ENV["GIT_SHA"]
    repository_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    try
        return readchomp(`git -C $repository_root rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function source_fingerprint()
    repository_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    paths = [
        joinpath(@__DIR__, "EffectorBSDE.jl"),
        joinpath(@__DIR__, "validate_sde_pipeline.jl"),
        joinpath(@__DIR__, "run_near_fold_pilot.jl"),
        joinpath(repository_root, "Project.toml"),
        joinpath(repository_root, "Manifest.toml"),
    ]
    payload = UInt8[]
    for path in paths
        append!(payload, codeunits(relpath(path, repository_root)))
        push!(payload, 0x00)
        append!(payload, read(path))
        push!(payload, 0x00)
    end
    return bytes2hex(sha256(payload))
end

"""
Numerical and observational settings for one stochastic trajectory.

`noise_mode = :paper_common` uses one shared Wiener process, matching the
explicit scalar `WienerProcess` in the original CIR notebook. The optional
`:independent` mode gives each population an independent Wiener process and
must be treated as a different stochastic model.

`gompertz_mode = :continuation` uses `a1*T*log(h/T)`, matching both the
original paper's stochastic ODE and the vector field used to calculate
`DEFAULT_B6_STAR`. The optional `:legacy_guarded` mode preserves the later
guard used by this repository's preliminary SDE scripts.
"""
Base.@kwdef struct SimulationConfig
    horizon::Float64 = 365.0
    saveat::Float64 = 1.0
    dt::Float64 = 0.01
    dtmax::Float64 = 0.05
    adaptive::Bool = true
    # SOSRI's strong-error controller is intentionally kept near the SciML
    # SDE defaults. ODE-style tolerances (for example 1e-5/1e-7) force
    # impractically small stochastic steps without improving this binary
    # establishment endpoint.
    reltol::Float64 = 1.0e-2
    abstol::Float64 = 1.0e-4
    maxiters::Int = 10_000_000
    noise_scale::Float64 = 1.0
    noise_mode::Symbol = :paper_common
    gompertz_mode::Symbol = :continuation
    paper_threshold::Float64 = 1.0
    strict_threshold::Float64 = 1.0e3
    strict_window_days::Float64 = 30.0
    strict_required_fraction::Float64 = 0.90
end

function validate_config(config::SimulationConfig)
    config.horizon > 0 || throw(ArgumentError("horizon must be positive"))
    config.saveat > 0 || throw(ArgumentError("saveat must be positive"))
    config.dt > 0 || throw(ArgumentError("dt must be positive"))
    config.dtmax >= config.dt ||
        throw(ArgumentError("dtmax must be greater than or equal to dt"))
    config.noise_scale >= 0 ||
        throw(ArgumentError("noise_scale must be nonnegative"))
    config.maxiters > 0 || throw(ArgumentError("maxiters must be positive"))
    config.noise_mode in (:paper_common, :independent) ||
        throw(ArgumentError("noise_mode must be :paper_common or :independent"))
    config.gompertz_mode in (:continuation, :legacy_guarded) ||
        throw(ArgumentError(
            "gompertz_mode must be :continuation or :legacy_guarded",
        ))
    0.0 <= config.strict_required_fraction <= 1.0 ||
        throw(ArgumentError("strict_required_fraction must be in [0, 1]"))
    return config
end

function base_parameters(;
    a10 = 2.0,
    b5 = 1.0e-4,
    a11 = 1.0e-1,
    b6 = DEFAULT_B6_STAR,
    b1 = 3.5e-6,
    b2 = 1.1e-7,
    noise_scale = 1.0,
    gompertz_mode = :continuation,
)
    return (
        a1 = 1.0e-1,
        h = 1.0e7,
        b1 = Float64(b1),
        b2 = Float64(b2),
        z1 = 0.0,
        a2 = 1.0e2,
        a3 = 1.0e8,
        z2 = 2.0e-1,
        a4 = 1.4e4,
        a5 = 2.5e-2,
        b3 = 4.0e-5,
        z3 = 4.12e-2,
        a6 = 1.1e-7,
        a7 = 1.0e-1,
        b4 = 1.0e-4,
        z4 = 2.0e-2,
        g1 = 1.0e10,
        g2 = 2.02e7,
        g3 = 2.02e7,
        a8 = 1.0e4,
        a9 = 5.0e-2,
        b5 = Float64(b5),
        z5 = 2.0e-2,
        g4 = 2.02e7,
        g5 = 1.0e3,
        a10 = Float64(a10),
        a11 = Float64(a11),
        b6 = Float64(b6),
        noise_scale = Float64(noise_scale),
        gompertz_mode = gompertz_mode,
    )
end

"""
The inoculation state used in the existing repository and in the original
paper's stochastic establishment experiment: two tumor cells, tumor-free
MDSC/NK baselines, zero CTL, and the corresponding B-cell baseline.

With direct B-to-CTL stimulation, CTL(0)=0 is an inoculation convention rather
than the exact five-state tumor-free equilibrium.
"""
function inoculation_start(p; tumor_cells = 2.0)
    xM0 = p.a2 / p.z2
    xN0 = p.z2 * p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC0 = 0.0
    xB0 = p.a8 / (p.z5 + p.b5 * xM0)
    return Float64[tumor_cells, xM0, xN0, xC0, xB0]
end

function gompertz_growth(xT, p)
    xT <= 0.0 && return 0.0
    if p.gompertz_mode === :continuation
        return p.a1 * xT * log(p.h / xT)
    elseif p.gompertz_mode === :legacy_guarded
        return p.a1 * xT * log(max(p.h / xT, 1.0))
    end
    throw(ArgumentError("unknown gompertz_mode $(p.gompertz_mode)"))
end

"""
Five-state deterministic drift used by both ODE and SDE problems.

The `max` guards only define trial evaluations outside the biological domain.
Accepted solution states are constrained separately through `isoutofdomain`;
the guards do not project accepted trajectories.
"""
function model_B!(du, u, p, t)
    xT = max(u[1], 0.0)
    xM = max(u[2], 0.0)
    xN = max(u[3], 0.0)
    xC = max(u[4], 0.0)
    xB = max(u[5], 0.0)

    b_to_ctl = p.a11 * xB / (p.g5 + xB)

    du[1] =
        gompertz_growth(xT, p) - p.b1 * xT * xN - p.b2 * xT * xC -
        p.b6 * xT * xB - p.z1 * xT
    du[2] = p.a2 + p.a3 * xT / (p.g1 + xT) - p.z2 * xM
    du[3] =
        p.a4 +
        p.a5 * xT^2 / (p.g2 + xT^2) *
        (1.0 + p.a10 * xB / (p.g5 + xB)) -
        p.b3 * xM * xN - p.z3 * xN
    du[4] =
        p.a6 * xT * xN + p.a7 * xT^2 / (p.g3 + xT^2) + b_to_ctl -
        p.b4 * xM * xC - p.z4 * xC
    du[5] =
        p.a8 + p.a9 * xT^2 / (p.g4 + xT^2) - p.b5 * xM * xB -
        p.z5 * xB

    # Define restorative drift during rejected trial evaluations outside the
    # positive domain, as recommended for nonnegative population models.
    @inbounds for i in eachindex(u)
        if u[i] < 0.0
            du[i] = max(du[i], 0.0)
        end
    end
    return nothing
end

function diffusion_common!(du, u, p, t)
    @inbounds for i in eachindex(u)
        du[i] = p.noise_scale * max(u[i], 0.0)
    end
    return nothing
end

function diffusion_independent!(du, u, p, t)
    @inbounds for i in eachindex(u)
        du[i] = p.noise_scale * max(u[i], 0.0)
    end
    return nothing
end

mutable struct EventTracker
    removed::Bool
    removal_time::Float64
end

function paper_removal_callback(tracker::EventTracker, threshold::Float64)
    condition(u, t, integrator) = u[1] < threshold
    function affect!(integrator)
        tracker.removed = true
        tracker.removal_time = Float64(integrator.t)
        terminate!(integrator)
    end
    return DiscreteCallback(
        condition,
        affect!;
        save_positions = (true, true),
    )
end

nonnegative_domain(u, p, t) = any(x -> !isfinite(x) || x < 0.0, u)

"""
Solve one SDE trajectory with an explicit per-trajectory seed.

A trajectory is permitted to stop before the horizon only through the tumor
removal event (`T < paper_threshold`). Any other incomplete or failed solve is
marked invalid later by `trajectory_metrics`.
"""
function solve_sde_once(;
    a10 = 2.0,
    b5 = 1.0e-4,
    a11 = 1.0e-1,
    b6 = DEFAULT_B6_STAR,
    b1 = 3.5e-6,
    b2 = 1.1e-7,
    tumor_cells = 2.0,
    seed::Integer = 1,
    config::SimulationConfig = SimulationConfig(),
)
    validate_config(config)
    p = base_parameters(
        a10 = a10,
        b5 = b5,
        a11 = a11,
        b6 = b6,
        b1 = b1,
        b2 = b2,
        noise_scale = config.noise_scale,
        gompertz_mode = config.gompertz_mode,
    )
    u0 = inoculation_start(p; tumor_cells = tumor_cells)
    tracker = EventTracker(false, NaN)
    callback = paper_removal_callback(tracker, config.paper_threshold)
    tspan = (0.0, config.horizon)

    if config.noise_mode === :paper_common
        noise = WienerProcess(0.0, 0.0, 0.0)
        prob = SDEProblem(
            model_B!,
            diffusion_common!,
            u0,
            tspan,
            p;
            noise = noise,
        )
    else
        prob = SDEProblem(model_B!, diffusion_independent!, u0, tspan, p)
    end

    sol = solve(
        prob,
        SOSRI();
        seed = Int(seed),
        callback = callback,
        dt = config.dt,
        dtmax = config.dtmax,
        adaptive = config.adaptive,
        saveat = config.saveat,
        reltol = config.reltol,
        abstol = config.abstol,
        maxiters = config.maxiters,
        isoutofdomain = nonnegative_domain,
    )

    return (
        solution = sol,
        removed = tracker.removed,
        removal_time = tracker.removal_time,
        seed = Int(seed),
        parameters = p,
        config = config,
    )
end

function solve_ode_reference(;
    a10 = 2.0,
    b5 = 1.0e-4,
    a11 = 1.0e-1,
    b6 = DEFAULT_B6_STAR,
    b1 = 3.5e-6,
    b2 = 1.1e-7,
    tumor_cells = 2.0,
    config::SimulationConfig = SimulationConfig(noise_scale = 0.0),
)
    validate_config(config)
    p = base_parameters(
        a10 = a10,
        b5 = b5,
        a11 = a11,
        b6 = b6,
        b1 = b1,
        b2 = b2,
        noise_scale = 0.0,
        gompertz_mode = config.gompertz_mode,
    )
    u0 = inoculation_start(p; tumor_cells = tumor_cells)
    prob = ODEProblem(model_B!, u0, (0.0, config.horizon), p)
    return solve(
        prob,
        Tsit5();
        saveat = config.saveat,
        reltol = min(config.reltol, 1.0e-9),
        abstol = min(config.abstol, 1.0e-11),
        isoutofdomain = nonnegative_domain,
    )
end

"""
Classify a saved tumor path under two observational definitions.

The paper-primary definition requires survival for the full horizon without a
recorded `T < 1` event. The stricter sensitivity definition additionally
requires tumor burden to stay above `strict_threshold` for at least the
configured fraction of the final time window.
"""
function classify_tumor_path(
    tumor,
    times;
    config::SimulationConfig = SimulationConfig(),
    removed_event::Bool = false,
)
    isempty(tumor) && throw(ArgumentError("tumor path cannot be empty"))
    length(tumor) == length(times) ||
        throw(ArgumentError("tumor and time arrays must have equal length"))

    horizon_tol = max(1.0e-8, 10eps(config.horizon))
    reached_horizon = times[end] >= config.horizon - horizon_tol
    crossed_from_saved_path = minimum(tumor) < config.paper_threshold
    paper_established =
        reached_horizon && !removed_event && !crossed_from_saved_path

    late_start = max(0.0, config.horizon - config.strict_window_days)
    late_idx = findall(t -> t >= late_start, times)
    late_fraction =
        isempty(late_idx) ? 0.0 :
        count(i -> tumor[i] >= config.strict_threshold, late_idx) /
        length(late_idx)
    strict_established =
        paper_established &&
        late_fraction >= config.strict_required_fraction

    strict_state =
        removed_event || crossed_from_saved_path ? "removed" :
        strict_established ? "strict_established" :
        paper_established ? "persistent_intermediate" : "invalid"

    return (
        reached_horizon = reached_horizon,
        paper_established = paper_established,
        paper_removed = removed_event || crossed_from_saved_path,
        strict_established = strict_established,
        strict_state = strict_state,
        late_fraction_above_strict = late_fraction,
    )
end

function first_saved_time_ge(values, times, threshold)
    idx = findfirst(x -> x >= threshold, values)
    return isnothing(idx) ? NaN : Float64(times[idx])
end

function solver_rejections(sol)
    if hasproperty(sol, :stats) && hasproperty(sol.stats, :nreject)
        return Int(sol.stats.nreject)
    elseif hasproperty(sol, :destats) && hasproperty(sol.destats, :nreject)
        return Int(sol.destats.nreject)
    end
    return -1
end

"""
Return auditable per-trajectory metrics.

Only a successful full-horizon solve or an explicit paper-removal termination
is valid. Solver failures and unexplained early endings remain in the raw CSV
but must not enter biological probability denominators.
"""
function trajectory_metrics(run)
    sol = run.solution
    config = run.config
    tumor = collect(sol[1, :])
    times = collect(sol.t)
    states = reduce(hcat, sol.u)
    classification = classify_tumor_path(
        tumor,
        times;
        config = config,
        removed_event = run.removed,
    )

    finite_states = all(isfinite, states)
    min_state = minimum(states)
    negativity_violation = min_state < 0.0
    retcode_text = string(sol.retcode)
    retcode_ok =
        successful_retcode(sol) || occursin("Terminated", retcode_text)
    explained_end = classification.reached_horizon || run.removed
    valid = retcode_ok && explained_end && finite_states && !negativity_violation

    return (
        seed = run.seed,
        retcode = retcode_text,
        valid = valid,
        termination_reason =
            run.removed ? "paper_removal" :
            classification.reached_horizon ? "horizon" : "solver_failure",
        t_end = Float64(times[end]),
        removal_time = run.removal_time,
        paper_established = valid && classification.paper_established,
        paper_removed = valid && classification.paper_removed,
        strict_established = valid && classification.strict_established,
        strict_state = valid ? classification.strict_state : "invalid",
        late_fraction_above_strict =
            classification.late_fraction_above_strict,
        T_min_saved = minimum(tumor),
        T_max_saved = maximum(tumor),
        T_final = Float64(states[1, end]),
        M_final = Float64(states[2, end]),
        NK_final = Float64(states[3, end]),
        CTL_final = Float64(states[4, end]),
        B_final = Float64(states[5, end]),
        saved_time_to_1e2 = first_saved_time_ge(tumor, times, 1.0e2),
        saved_time_to_1e3 = first_saved_time_ge(tumor, times, 1.0e3),
        saved_time_to_1e5 = first_saved_time_ge(tumor, times, 1.0e5),
        min_state_saved = Float64(min_state),
        negativity_violation = negativity_violation,
        solver_rejections = solver_rejections(sol),
    )
end

end
