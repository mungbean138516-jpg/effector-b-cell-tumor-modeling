using Test

include(
    joinpath(
        @__DIR__,
        "..",
        "code",
        "post_meeting",
        "sde_validation",
        "EffectorBSDE.jl",
    ),
)
using .EffectorBSDE

@testset "SDE criterion classifiers" begin
    config = SimulationConfig(
        horizon = 365.0,
        strict_threshold = 1.0e3,
        strict_window_days = 30.0,
        strict_required_fraction = 0.90,
    )
    times = collect(0.0:1.0:365.0)

    survived = fill(2.0, length(times))
    survived[times .>= 335.0] .= 2.0e3
    result = classify_tumor_path(survived, times; config = config)
    @test result.paper_established
    @test result.strict_established
    @test result.strict_state == "strict_established"

    dipped_then_recovered = fill(2.0e3, length(times))
    dipped_then_recovered[11] = 0.9
    result =
        classify_tumor_path(dipped_then_recovered, times; config = config)
    @test !result.paper_established
    @test result.paper_removed
    @test result.strict_state == "removed"

    persistent_small = fill(2.0, length(times))
    result = classify_tumor_path(persistent_small, times; config = config)
    @test result.paper_established
    @test !result.strict_established
    @test result.strict_state == "persistent_intermediate"

    incomplete_times = collect(0.0:1.0:100.0)
    incomplete = fill(2.0e3, length(incomplete_times))
    result =
        classify_tumor_path(incomplete, incomplete_times; config = config)
    @test !result.reached_horizon
    @test !result.paper_established
    @test result.strict_state == "invalid"
end

@testset "Model and solver smoke checks" begin
    config = SimulationConfig(
        horizon = 10.0,
        saveat = 0.5,
        dt = 0.01,
        dtmax = 0.05,
        noise_scale = 0.0,
    )
    p = base_parameters(
        b6 = 0.8 * DEFAULT_B6_STAR,
        noise_scale = 0.0,
    )
    u0 = inoculation_start(p)
    du = similar(u0)
    model_B!(du, u0, p, 0.0)
    @test all(isfinite, du)
    @test all(u0 .>= 0)

    run_1 = solve_sde_once(
        b6 = 0.8 * DEFAULT_B6_STAR,
        seed = 11,
        config = config,
    )
    run_2 = solve_sde_once(
        b6 = 0.8 * DEFAULT_B6_STAR,
        seed = 12,
        config = config,
    )
    metrics = trajectory_metrics(run_1)
    @test metrics.valid
    @test !metrics.negativity_violation
    @test run_1.solution.t == run_2.solution.t
    @test run_1.solution.u == run_2.solution.u

    independent_run = solve_sde_once(
        seed = 13,
        config = SimulationConfig(
            horizon = 2.0,
            saveat = 0.5,
            noise_scale = 0.1,
            noise_mode = :independent,
        ),
    )
    @test trajectory_metrics(independent_run).valid

    removal_run = solve_sde_once(
        b6 = 0.8 * DEFAULT_B6_STAR,
        seed = 11,
        config = SimulationConfig(
            horizon = 5.0,
            saveat = 0.5,
            noise_scale = 1.0,
            noise_mode = :paper_common,
        ),
    )
    removal_metrics = trajectory_metrics(removal_run)
    @test removal_metrics.valid
    @test removal_run.removed
    @test removal_metrics.termination_reason == "paper_removal"
    @test !removal_metrics.paper_established

    interrupted_run = solve_sde_once(
        seed = 17,
        config = SimulationConfig(
            horizon = 10.0,
            saveat = 0.5,
            noise_scale = 0.1,
            maxiters = 1,
        ),
    )
    interrupted_metrics = trajectory_metrics(interrupted_run)
    @test !interrupted_metrics.valid
    @test interrupted_metrics.termination_reason == "solver_failure"
    @test !interrupted_metrics.paper_established
end
