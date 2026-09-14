using Test

include(
    joinpath(
        @__DIR__,
        "..",
        "code",
        "post_meeting",
        "attractor_switching_poc",
        "AttractorSwitchingPOC.jl",
    ),
)
using .AttractorSwitchingPOC

const CONTINUATION_PATH = joinpath(
    @__DIR__,
    "..",
    "results",
    "post_meeting",
    "formal_b6_continuation",
    "formal_continuation_all_points.csv",
)

@testset "Attractor geometry" begin
    attractors = refine_attractor_set(CONTINUATION_PATH; b6 = 1.5e-6)
    @test attractors.low[1] < attractors.saddle[1] < attractors.high[1]
    @test attractors.low_n_unstable == 0
    @test attractors.saddle_n_unstable == 1
    @test attractors.high_n_unstable == 0
    @test maximum([
        attractors.low_residual,
        attractors.saddle_residual,
        attractors.high_residual,
    ]) < 1.0e-9

    low_class = classify_relaxed_state(attractors.low, attractors)
    high_class = classify_relaxed_state(attractors.high, attractors)
    @test low_class.classification == "low"
    @test high_class.classification == "high"
    @test log_state_distance(attractors.low, attractors.low) < 1.0e-12
end

@testset "Switching solver invariants" begin
    attractors = refine_attractor_set(CONTINUATION_PATH; b6 = 1.5e-6)
    config = SwitchingConfig(
        noise_horizon = 5.0,
        relax_horizon = 200.0,
        saveat = 0.5,
        relax_saveat = 1.0,
    )

    low_zero = run_switching_once(
        :low,
        attractors;
        noise_scale = 0.0,
        seed = 101,
        config = config,
    )
    high_zero = run_switching_once(
        :high,
        attractors;
        noise_scale = 0.0,
        seed = 102,
        config = config,
    )
    @test low_zero.valid
    @test high_zero.valid
    @test !low_zero.switched
    @test !high_zero.switched
    @test low_zero.final_attractor == "low"
    @test high_zero.final_attractor == "high"

    stochastic_1 = solve_noise_pulse(
        attractors.low;
        noise_scale = 0.10,
        seed = 303,
        config = config,
    )
    stochastic_2 = solve_noise_pulse(
        attractors.low;
        noise_scale = 0.10,
        seed = 303,
        config = config,
    )
    @test stochastic_1.t == stochastic_2.t
    @test stochastic_1.u == stochastic_2.u
    @test minimum(reduce(hcat, stochastic_1.u)) >= 0.0

    metrics = run_to_metrics(low_zero, config)
    @test metrics.valid
    @test metrics.noise_reached_horizon
    @test metrics.relaxation_reached_horizon
end
