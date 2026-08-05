module AttractorSwitchingPOC

using CSV
using DataFrames
using ForwardDiff
using LinearAlgebra
using OrdinaryDiffEq: Tsit5
using SHA: sha256
using SciMLBase: ODEProblem, SDEProblem, solve, successful_retcode
using StochasticDiffEq: SOSRI, WienerProcess

include(joinpath(@__DIR__, "..", "sde_validation", "EffectorBSDE.jl"))
using .EffectorBSDE

export AttractorSet,
       SwitchingConfig,
       classify_relaxed_state,
       continuation_guesses,
       equilibrium_table,
       log_state_distance,
       poc_source_fingerprint,
       refine_attractor_set,
       repository_git_sha,
       run_switching_once,
       run_to_metrics,
       solve_noise_pulse

const STATE_NAMES = (:T, :M, :NK, :CTL, :B)
const STATE_SCALE = Float64[1.0e7, 5.0e5, 2.0e5, 1.0e3, 1.0e5]

Base.@kwdef struct SwitchingConfig
    b6::Float64 = 1.5e-6
    noise_horizon::Float64 = 60.0
    relax_horizon::Float64 = 1_500.0
    saveat::Float64 = 0.5
    relax_saveat::Float64 = 1.0
    dt::Float64 = 0.01
    dtmax::Float64 = 0.05
    adaptive::Bool = true
    reltol::Float64 = 1.0e-2
    abstol::Float64 = 1.0e-4
    maxiters::Int = 10_000_000
    classification_tolerance::Float64 = 0.03
end

struct AttractorSet
    low::Vector{Float64}
    saddle::Vector{Float64}
    high::Vector{Float64}
    low_residual::Float64
    saddle_residual::Float64
    high_residual::Float64
    low_max_real_eig::Float64
    saddle_max_real_eig::Float64
    high_max_real_eig::Float64
    low_n_unstable::Int
    saddle_n_unstable::Int
    high_n_unstable::Int
end

repository_root() = normpath(joinpath(@__DIR__, "..", "..", ".."))
repository_git_sha() = EffectorBSDE.repository_git_sha()

function poc_source_fingerprint()
    root = repository_root()
    paths = [
        joinpath(@__DIR__, "AttractorSwitchingPOC.jl"),
        joinpath(@__DIR__, "run_attractor_switching_poc.jl"),
        joinpath(@__DIR__, "..", "sde_validation", "EffectorBSDE.jl"),
        joinpath(root, "Project.toml"),
        joinpath(root, "Manifest.toml"),
        joinpath(
            root,
            "results",
            "post_meeting",
            "formal_b6_continuation",
            "formal_continuation_all_points.csv",
        ),
    ]
    payload = UInt8[]
    for path in paths
        isfile(path) || error("Fingerprint input is missing: $(path)")
        append!(payload, codeunits(relpath(path, root)))
        push!(payload, 0x00)
        append!(payload, read(path))
        push!(payload, 0x00)
    end
    return bytes2hex(sha256(payload))
end

function validate_config(config::SwitchingConfig)
    config.b6 > 0.0 || throw(ArgumentError("b6 must be positive"))
    config.noise_horizon > 0.0 ||
        throw(ArgumentError("noise_horizon must be positive"))
    config.relax_horizon > 0.0 ||
        throw(ArgumentError("relax_horizon must be positive"))
    config.saveat > 0.0 || throw(ArgumentError("saveat must be positive"))
    config.relax_saveat > 0.0 ||
        throw(ArgumentError("relax_saveat must be positive"))
    config.dt > 0.0 || throw(ArgumentError("dt must be positive"))
    config.dtmax >= config.dt ||
        throw(ArgumentError("dtmax must be at least dt"))
    config.maxiters > 0 || throw(ArgumentError("maxiters must be positive"))
    config.classification_tolerance > 0.0 ||
        throw(ArgumentError("classification_tolerance must be positive"))
    return config
end

function state_from_row(row)
    return Float64[row.T, row.M, row.NK, row.CTL, row.B]
end

function nearest_row(frame::DataFrame, target_q::Float64, label::AbstractString)
    nrow(frame) > 0 || error("No continuation rows found for $(label)")
    index = argmin(abs.(Float64.(frame.q) .- target_q))
    return frame[index, :]
end

"""
Select low-attractor, saddle, and high-attractor starting guesses from the
formal continuation export. The guesses are refined at the exact requested
`b6` before any stochastic simulation is run.
"""
function continuation_guesses(path::AbstractString, b6::Real)
    data = CSV.read(path, DataFrame)
    required = Set([:seed, :q, :T, :M, :NK, :CTL, :B, :stable])
    missing_columns = setdiff(required, Set(propertynames(data)))
    isempty(missing_columns) ||
        error("Continuation table is missing columns: $(missing_columns)")

    target_q = log10(Float64(b6))
    high_rows = data[(data.seed .== "high_seed") .& data.stable, :]
    low_rows = data[(data.seed .== "low_seed") .& data.stable, :]
    saddle_rows = data[(data.seed .== "low_seed") .& .!data.stable, :]

    high_row = nearest_row(high_rows, target_q, "high attractor")
    low_row = nearest_row(low_rows, target_q, "low attractor")
    saddle_row = nearest_row(saddle_rows, target_q, "saddle")
    guesses = (
        low = state_from_row(low_row),
        saddle = state_from_row(saddle_row),
        high = state_from_row(high_row),
    )
    guesses.low[1] < guesses.saddle[1] < guesses.high[1] || error(
        "Continuation guesses are not ordered low < saddle < high in tumor burden",
    )
    return guesses
end

function drift_vector(x::AbstractVector, p)
    dx = similar(x)
    EffectorBSDE.model_B!(dx, x, p, 0.0)
    return dx
end

state_from_log(y::AbstractVector) = STATE_SCALE .* exp.(y)
log_from_state(x::AbstractVector) = log.(x ./ STATE_SCALE)

function scaled_equilibrium_residual(y::AbstractVector, p)
    x = state_from_log(y)
    return drift_vector(x, p) ./ x
end

function refine_equilibrium(
    initial_state::AbstractVector,
    p;
    tolerance::Float64 = 1.0e-11,
    max_iterations::Int = 50,
)
    all(initial_state .> 0.0) ||
        throw(ArgumentError("equilibrium guess must be strictly positive"))
    y = log_from_state(Float64.(initial_state))
    best_y = copy(y)
    best_residual = norm(scaled_equilibrium_residual(y, p), Inf)
    converged = best_residual < tolerance
    iterations = 0

    for iteration in 1:max_iterations
        iterations = iteration
        residual = scaled_equilibrium_residual(y, p)
        residual_norm = norm(residual, Inf)
        if residual_norm < best_residual
            best_y = copy(y)
            best_residual = residual_norm
        end
        if residual_norm < tolerance
            converged = true
            break
        end

        jacobian = ForwardDiff.jacobian(
            yy -> scaled_equilibrium_residual(yy, p),
            y,
        )
        step = try
            -(jacobian \ residual)
        catch
            -((jacobian + 1.0e-10I) \ residual)
        end
        largest_step = maximum(abs.(step))
        largest_step > 1.0 && (step .*= 1.0 / largest_step)

        accepted = false
        alpha = 1.0
        for _ in 1:20
            candidate = y + alpha * step
            candidate_norm = norm(
                scaled_equilibrium_residual(candidate, p),
                Inf,
            )
            if isfinite(candidate_norm) && candidate_norm < residual_norm
                y = candidate
                accepted = true
                break
            end
            alpha *= 0.5
        end
        accepted || break
    end

    final_state = state_from_log(best_y)
    eigenvalues = eigvals(
        ForwardDiff.jacobian(xx -> drift_vector(xx, p), final_state),
    )
    max_real_eig = maximum(real.(eigenvalues))
    n_unstable = count(real.(eigenvalues) .> 1.0e-8)
    return (
        state = Float64.(final_state),
        residual = Float64(best_residual),
        converged = converged || best_residual < tolerance,
        iterations = iterations,
        eigenvalues = eigenvalues,
        max_real_eig = Float64(max_real_eig),
        n_unstable = n_unstable,
    )
end

"""Refine and verify the three equilibria used by the switching experiment."""
function refine_attractor_set(
    continuation_path::AbstractString;
    b6::Float64 = 1.5e-6,
)
    guesses = continuation_guesses(continuation_path, b6)
    p = EffectorBSDE.base_parameters(
        b6 = b6,
        noise_scale = 0.0,
        gompertz_mode = :continuation,
    )
    low = refine_equilibrium(guesses.low, p)
    saddle = refine_equilibrium(guesses.saddle, p)
    high = refine_equilibrium(guesses.high, p)

    all((low.converged, saddle.converged, high.converged)) ||
        error("At least one equilibrium refinement failed")
    low.state[1] < saddle.state[1] < high.state[1] ||
        error("Refined equilibria are not ordered low < saddle < high")
    low.n_unstable == 0 || error("Low-tumor equilibrium is not stable")
    saddle.n_unstable == 1 || error("Middle equilibrium is not a saddle")
    high.n_unstable == 0 || error("High-tumor equilibrium is not stable")

    return AttractorSet(
        low.state,
        saddle.state,
        high.state,
        low.residual,
        saddle.residual,
        high.residual,
        low.max_real_eig,
        saddle.max_real_eig,
        high.max_real_eig,
        low.n_unstable,
        saddle.n_unstable,
        high.n_unstable,
    )
end

function equilibrium_table(attractors::AttractorSet, b6::Real)
    rows = DataFrame()
    for label in (:low, :saddle, :high)
        state = getproperty(attractors, label)
        push!(
            rows,
            (
                equilibrium = string(label),
                b6 = Float64(b6),
                T = state[1],
                M = state[2],
                NK = state[3],
                CTL = state[4],
                B = state[5],
                residual = getproperty(attractors, Symbol(label, :_residual)),
                max_real_eig =
                    getproperty(attractors, Symbol(label, :_max_real_eig)),
                n_unstable =
                    getproperty(attractors, Symbol(label, :_n_unstable)),
                stable = getproperty(
                    attractors,
                    Symbol(label, :_n_unstable),
                ) == 0,
            ),
        )
    end
    return rows
end

function log_state_distance(x::AbstractVector, target::AbstractVector)
    length(x) == length(STATE_SCALE) ||
        throw(ArgumentError("state has the wrong dimension"))
    length(target) == length(STATE_SCALE) ||
        throw(ArgumentError("target has the wrong dimension"))
    floor_values = 1.0e-12 .* STATE_SCALE
    safe_x = max.(Float64.(x), floor_values)
    safe_target = max.(Float64.(target), floor_values)
    return norm(log10.(safe_x ./ safe_target)) / sqrt(length(STATE_SCALE))
end

function classify_relaxed_state(
    state::AbstractVector,
    attractors::AttractorSet;
    tolerance::Float64 = 0.03,
)
    low_distance = log_state_distance(state, attractors.low)
    high_distance = log_state_distance(state, attractors.high)
    nearest = low_distance <= high_distance ? "low" : "high"
    nearest_distance = min(low_distance, high_distance)
    classification = nearest_distance <= tolerance ? nearest : "ambiguous"
    return (
        classification = classification,
        low_distance = low_distance,
        high_distance = high_distance,
        nearest_distance = nearest_distance,
    )
end

function solve_noise_pulse(
    initial_state::AbstractVector;
    noise_scale::Float64,
    seed::Integer,
    config::SwitchingConfig = SwitchingConfig(),
)
    validate_config(config)
    noise_scale >= 0.0 ||
        throw(ArgumentError("noise_scale must be nonnegative"))
    p = EffectorBSDE.base_parameters(
        b6 = config.b6,
        noise_scale = noise_scale,
        gompertz_mode = :continuation,
    )
    noise = WienerProcess(0.0, 0.0, 0.0)
    problem = SDEProblem(
        EffectorBSDE.model_B!,
        EffectorBSDE.diffusion_common!,
        Float64.(initial_state),
        (0.0, config.noise_horizon),
        p;
        noise = noise,
    )
    return solve(
        problem,
        SOSRI();
        seed = Int(seed),
        dt = config.dt,
        dtmax = config.dtmax,
        adaptive = config.adaptive,
        saveat = config.saveat,
        reltol = config.reltol,
        abstol = config.abstol,
        maxiters = config.maxiters,
        isoutofdomain = EffectorBSDE.nonnegative_domain,
    )
end

function solve_deterministic_relaxation(
    initial_state::AbstractVector,
    config::SwitchingConfig,
)
    p = EffectorBSDE.base_parameters(
        b6 = config.b6,
        noise_scale = 0.0,
        gompertz_mode = :continuation,
    )
    problem = ODEProblem(
        EffectorBSDE.model_B!,
        Float64.(initial_state),
        (0.0, config.relax_horizon),
        p,
    )
    return solve(
        problem,
        Tsit5();
        saveat = config.relax_saveat,
        reltol = 1.0e-9,
        abstol = 1.0e-11,
        maxiters = config.maxiters,
        isoutofdomain = EffectorBSDE.nonnegative_domain,
    )
end

function solution_diagnostics(solution, expected_horizon::Float64)
    states = reduce(hcat, solution.u)
    finite = all(isfinite, states)
    minimum_state = finite ? minimum(states) : -Inf
    nonnegative = minimum_state >= 0.0
    reached_horizon =
        solution.t[end] >= expected_horizon - max(1.0e-8, 10eps(expected_horizon))
    valid =
        successful_retcode(solution) && reached_horizon && finite && nonnegative
    return (
        valid = valid,
        retcode = string(solution.retcode),
        reached_horizon = reached_horizon,
        finite = finite,
        nonnegative = nonnegative,
        minimum_state = Float64(minimum_state),
        solver_rejections = hasproperty(solution.stats, :nreject) ?
                            Int(solution.stats.nreject) : -1,
    )
end

function run_switching_once(
    initial_attractor::Symbol,
    attractors::AttractorSet;
    noise_scale::Float64,
    seed::Integer,
    config::SwitchingConfig = SwitchingConfig(),
)
    initial_attractor in (:low, :high) ||
        throw(ArgumentError("initial_attractor must be :low or :high"))
    initial_state = getproperty(attractors, initial_attractor)
    noise_solution = solve_noise_pulse(
        initial_state;
        noise_scale = noise_scale,
        seed = seed,
        config = config,
    )
    noise_diagnostics =
        solution_diagnostics(noise_solution, config.noise_horizon)

    if !noise_diagnostics.valid
        return (
            initial_attractor = string(initial_attractor),
            final_attractor = "invalid",
            switched = false,
            valid = false,
            seed = Int(seed),
            noise_scale = noise_scale,
            noise_solution = noise_solution,
            relaxation_solution = nothing,
            noise_diagnostics = noise_diagnostics,
            relaxation_diagnostics = nothing,
            classification = nothing,
        )
    end

    relaxation_solution =
        solve_deterministic_relaxation(noise_solution.u[end], config)
    relaxation_diagnostics =
        solution_diagnostics(relaxation_solution, config.relax_horizon)
    classification = classify_relaxed_state(
        relaxation_solution.u[end],
        attractors;
        tolerance = config.classification_tolerance,
    )
    valid =
        relaxation_diagnostics.valid &&
        classification.classification != "ambiguous"
    final_attractor = valid ? classification.classification : "ambiguous"
    switched = valid && final_attractor != string(initial_attractor)

    return (
        initial_attractor = string(initial_attractor),
        final_attractor = final_attractor,
        switched = switched,
        valid = valid,
        seed = Int(seed),
        noise_scale = noise_scale,
        noise_solution = noise_solution,
        relaxation_solution = relaxation_solution,
        noise_diagnostics = noise_diagnostics,
        relaxation_diagnostics = relaxation_diagnostics,
        classification = classification,
    )
end

function state_metrics(prefix::Symbol, state::AbstractVector)
    return NamedTuple{
        Tuple(Symbol(prefix, :_, name) for name in STATE_NAMES)
    }(Tuple(Float64.(state)))
end

function run_to_metrics(run, config::SwitchingConfig)
    noise_states = reduce(hcat, run.noise_solution.u)
    noise_endpoint = run.noise_solution.u[end]
    relaxation_endpoint =
        isnothing(run.relaxation_solution) ? fill(NaN, 5) :
        run.relaxation_solution.u[end]
    relaxation_diagnostics = run.relaxation_diagnostics
    classification = run.classification

    return merge(
        (
            initial_attractor = run.initial_attractor,
            final_attractor = run.final_attractor,
            switched = run.switched,
            valid = run.valid,
            seed = run.seed,
            b6 = config.b6,
            noise_scale = run.noise_scale,
            noise_horizon = config.noise_horizon,
            relax_horizon = config.relax_horizon,
            noise_retcode = run.noise_diagnostics.retcode,
            relaxation_retcode = isnothing(relaxation_diagnostics) ?
                                  "not_run" :
                                  relaxation_diagnostics.retcode,
            noise_reached_horizon = run.noise_diagnostics.reached_horizon,
            relaxation_reached_horizon =
                isnothing(relaxation_diagnostics) ? false :
                relaxation_diagnostics.reached_horizon,
            min_state_saved = minimum(noise_states),
            T_min_noise = minimum(noise_states[1, :]),
            T_max_noise = maximum(noise_states[1, :]),
            low_distance = isnothing(classification) ? NaN :
                           classification.low_distance,
            high_distance = isnothing(classification) ? NaN :
                            classification.high_distance,
            nearest_distance = isnothing(classification) ? NaN :
                               classification.nearest_distance,
            noise_solver_rejections =
                run.noise_diagnostics.solver_rejections,
            relaxation_solver_rejections =
                isnothing(relaxation_diagnostics) ? -1 :
                relaxation_diagnostics.solver_rejections,
        ),
        state_metrics(:noise_end, noise_endpoint),
        state_metrics(:relax_end, relaxation_endpoint),
    )
end

end
