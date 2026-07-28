# Compatibility entry point for the validated SDE workflow.
#
# The former implementation is preserved as
# `legacy_b6_sde_threshold_refined.jl` for provenance. New analyses must use
# the validation gates and raw per-trajectory QC in the post-meeting pipeline.

include(
    joinpath(
        @__DIR__,
        "..",
        "post_meeting",
        "sde_validation",
        "run_near_fold_pilot.jl",
    ),
)

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
