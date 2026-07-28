using BifurcationKit
using Accessors
using ForwardDiff
using LinearAlgebra
using CSV
using DataFrames
using Plots

const BK = BifurcationKit

# ------------------------------------------------------------
# Two-parameter continuation of the saddle-node in the b5-b6 plane
# ------------------------------------------------------------
# Scientific question:
#   How does MDSC suppression of effector B cells (b5) shift the
#   saddle-node threshold in B-cell anti-tumor strength (b6)?
#
# Method:
#   1. Reconstruct the one-parameter equilibrium branch in q6 = log10(b6)
#      at baseline q5 = log10(b5) = -4.
#   2. Select the detected zero-eigenvalue / fold point near the previously
#      classified saddle-node.
#   3. Continue that fold in the second parameter q5 = log10(b5) using
#      BifurcationKit's minimally augmented codimension-two continuation.
#   4. Export the formal fold curve and compare it with the earlier
#      simulation-defined T(10000) < 1e3 control boundary when available.
#
# Figure titles are descriptive and result-based only.
# ------------------------------------------------------------

outdir = "two_parameter_fold_continuation_outputs"
mkpath(outdir)

# ----------------------------
# User settings
# ----------------------------
Q5_BASE = -4.0                 # b5 = 1e-4
Q5_MIN  = -5.0                 # b5 = 1e-5
Q5_MAX  = -3.0                 # b5 = 1e-3

Q6_MIN  = -7.2                 # one-parameter seed branch range
Q6_MAX  = -5.2

LOW_SEED_B6 = 1.5e-6
FOLD_B6_REFERENCE = 6.676254e-7

# Existing state scaling used in the formal one-parameter continuation.
const SCALE = [1.0e7, 5.0e5, 2.0e5, 1.0e3, 1.0e5]

T_CONTROL = 1.0e3
T_DOMINANT = 1.0e5

# ----------------------------
# Locate existing inputs
# ----------------------------
function first_existing(paths; required=true)
    for p in paths
        if isfile(p)
            return p
        end
    end
    if required
        error("Could not find any expected input file:\n" * join(paths, "\n"))
    end
    return nothing
end

low_seed_file = first_existing([
    "formal_continuation_low_seed.csv",
    joinpath("formal_b6_continuation_outputs", "formal_continuation_low_seed.csv"),
    joinpath("results", "post_meeting", "formal_b6_continuation",
             "formal_continuation_low_seed.csv"),
    joinpath("results", "post_meeting", "formal_b6_continuation",
             "formal_b6_continuation_outputs", "formal_continuation_low_seed.csv")
])

classification_file = first_existing([
    "zero_eigenvalue_classification.csv",
    joinpath("zero_eigenvalue_classification_outputs",
             "zero_eigenvalue_classification.csv"),
    joinpath("results", "post_meeting", "zero_eigenvalue_classification",
             "zero_eigenvalue_classification.csv")
]; required=false)

empirical_boundary_file = first_existing([
    "step3_control_boundary.csv",
    joinpath("step3_b5_b6_outputs", "step3_control_boundary.csv"),
    joinpath("results", "post_meeting", "step3_heatmap", "raw",
             "step3_control_boundary.csv"),
    joinpath("results", "post_meeting", "step3_heatmap", "fixed",
             "step3_fixed_control_boundary.csv")
]; required=false)

low_seed_rows = CSV.read(low_seed_file, DataFrame)
println("Using low-equilibrium seed table: ", low_seed_file)

# Refined fold location from the previous local classification when available.
q6_fold_reference = log10(FOLD_B6_REFERENCE)

if classification_file !== nothing
    class_df = CSV.read(classification_file, DataFrame)
    lookup = Dict(string(r.quantity) => r.value for r in eachrow(class_df))
    if haskey(lookup, "b6_star")
        q6_fold_reference = log10(parse(Float64, string(lookup["b6_star"])))
    end
    println("Using refined saddle-node reference from: ", classification_file)
else
    println("Using fallback saddle-node reference b6 = ", FOLD_B6_REFERENCE)
end

# ----------------------------
# Model
# ----------------------------
function base_p(; q5=Q5_BASE, q6=-6.0)
    b5 = 10.0^q5
    b6 = 10.0^q6

    return (
        a1 = 1.0e-1,
        h  = 1.0e7,
        b1 = 3.5e-6,
        b2 = 1.1e-7,
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

        a8  = 1.0e4,
        a9  = 5.0e-2,
        b5  = b5,
        z5  = 2.0e-2,
        g4  = 2.02e7,
        g5  = 1.0e3,
        a10 = 2.0,
        a11 = 1.0e-1,
        b6  = b6
    )
end

function rhs_x(x, p)
    xT, xM, xN, xC, xB = x

    b_to_ctl = p.a11 * xB / (p.g5 + xB)
    growth = p.a1 * xT * log(p.h / xT)

    return [
        growth - p.b1*xT*xN - p.b2*xT*xC - p.b6*xT*xB - p.z1*xT,
        p.a2 + p.a3*xT/(p.g1 + xT) - p.z2*xM,
        p.a4 + p.a5*xT^2/(p.g2 + xT^2) *
            (1.0 + p.a10*xB/(p.g5 + xB)) -
            p.b3*xM*xN - p.z3*xN,
        p.a6*xT*xN + p.a7*xT^2/(p.g3 + xT^2) +
            b_to_ctl - p.b4*xM*xC - p.z4*xC,
        p.a8 + p.a9*xT^2/(p.g4 + xT^2) -
            p.b5*xM*xB - p.z5*xB
    ]
end

x_from_y(y) = SCALE .* exp.(y)
y_from_x(x) = log.(x ./ SCALE)

# Parameters passed to BifurcationKit are (q6, q5).
function F(y, par)
    x = x_from_y(y)
    p = base_p(q5=par.q5, q6=par.q6)
    return rhs_x(x, p) ./ x
end

function J(y, par)
    return ForwardDiff.jacobian(yy -> F(yy, par), y)
end

function record_equilibrium(y, par; k...)
    # In one-parameter continuation, BifurcationKit may pass the active
    # continuation parameter itself as `par` rather than the full parameter tuple.
    # Here the active parameter is q6; q5 stays fixed at Q5_BASE.
    q6 = par isa Real ? Float64(par) : Float64(par.q6)
    q5 = par isa Real ? Float64(Q5_BASE) : Float64(par.q5)

    x = x_from_y(y)
    return (
        q6 = q6,
        q5 = q5,
        b6 = 10.0^q6,
        b5 = 10.0^q5,
        T = x[1],
        M = x[2],
        NK = x[3],
        CTL = x[4],
        B = x[5],
        logT = log10(x[1] + 1.0)
    )
end

# ----------------------------
# Select a low-tumor equilibrium seed
# ----------------------------
idx_seed = argmin(abs.(log10.(low_seed_rows.b6) .- log10(LOW_SEED_B6)))
seed_row = low_seed_rows[idx_seed, :]

x_seed = Float64[
    seed_row.T,
    seed_row.M,
    seed_row.NK,
    seed_row.CTL,
    seed_row.B
]

y_seed = y_from_x(x_seed)
q6_seed = log10(Float64(seed_row.b6))
par_seed = (q6=q6_seed, q5=Q5_BASE)

println("\nSelected one-parameter continuation seed:")
println("  b5 = ", 10.0^Q5_BASE)
println("  b6 = ", 10.0^q6_seed)
println("  tumor = ", x_seed[1])

# ----------------------------
# One-parameter continuation in q6 at fixed q5
# ----------------------------
prob = BifurcationProblem(
    F,
    y_seed,
    par_seed,
    (@optic _.q6);
    J=J,
    record_from_solution=record_equilibrium
)

newton_par = NewtonPar(
    tol=1.0e-10,
    max_iterations=30,
    linesearch=true
)

one_parameter_options = ContinuationPar(
    ds=-0.005,
    dsmin=1.0e-5,
    dsmax=0.02,
    p_min=Q6_MIN,
    p_max=Q6_MAX,
    max_steps=1200,
    newton_options=newton_par,
    nev=5,
    save_eig_every_step=1,
    detect_bifurcation=3,
    detect_fold=true
)

println("\nReconstructing the one-parameter equilibrium branch...")
br = continuation(
    prob,
    PALC(tangent=Bordered()),
    one_parameter_options;
    bothside=true
)

open(joinpath(outdir, "one_parameter_branch_summary.txt"), "w") do io
    show(io, br)
end

println("\nDetected one-parameter special points:")
for (i, sp) in enumerate(br.specialpoint)
    println("  #", i,
            " type=", sp.type,
            " q6=", sp.param,
            " b6=", 10.0^sp.param)
end

# The equilibrium detector may label a fold as :bp; the codim-2 continuation
# routine treats any non-Hopf steady-state singularity as a fold candidate.
candidate_indices = [
    i for i in eachindex(br.specialpoint)
    if br.specialpoint[i].type != :endpoint
]

isempty(candidate_indices) &&
    error("No steady-state special point was detected on the one-parameter branch.")

indfold = candidate_indices[
    argmin(abs.([
        br.specialpoint[i].param - q6_fold_reference
        for i in candidate_indices
    ]))
]

sp0 = br.specialpoint[indfold]

println("\nSelected fold candidate:")
println("  special-point index = ", indfold)
println("  detector label      = ", sp0.type)
println("  q6                  = ", sp0.param)
println("  b6                  = ", 10.0^sp0.param)

# ----------------------------
# Continue the fold in q5 = log10(b5)
# ----------------------------
two_parameter_options = ContinuationPar(
    ds=0.01,
    dsmin=1.0e-5,
    dsmax=0.035,
    p_min=Q5_MIN,
    p_max=Q5_MAX,
    max_steps=1600,
    newton_options=newton_par,

    # Codimension-two special points are detected through events.
    detect_bifurcation=0,
    detect_event=2,
    nev=5
)

# Extract the original five-dimensional equilibrium state from the object
# passed by BifurcationKit's fold recorder.  Depending on the internal
# continuation problem, `z` may be either a BorderedArray with a `.u` field
# or already a 5-component view/vector of the equilibrium state.
function fold_state_y(z)
    nstate = length(SCALE)

    if hasproperty(z, :u)
        yraw = z.u
    elseif z isa AbstractVector
        if length(z) == nstate
            yraw = z
        elseif length(z) > nstate
            yraw = @view z[1:nstate]
        else
            error("record_fold received a state vector with length $(length(z)); expected at least $nstate.")
        end
    else
        error("record_fold received unsupported state object of type $(typeof(z)).")
    end

    return Float64.(collect(yraw))
end

# Record the full physical parameter/state values along the fold curve.
function record_fold(z, p2; iter, state, k...)
    pars = BK.getparams(iter, state)

    y = fold_state_y(z)
    x = x_from_y(y)

    Jstate = J(y, pars)
    S = svd(Jstate)
    v = S.V[:, end]
    w = S.U[:, end]

    if dot(w, v) < 0
        w .*= -1.0
    end

    v ./= norm(v)
    w ./= dot(w, v)

    hv = 1.0e-4
    Hvv = (
        F(y + hv*v, pars) -
        2.0F(y, pars) +
        F(y - hv*v, pars)
    ) / hv^2

    b20 = dot(w, Hvv)
    residual = norm(F(y, pars), Inf)

    eigs = eigvals(Jstate)
    eig_order = sortperm(abs.(eigs))
    nonzero_eigs = eigs[eig_order[2:end]]
    leading_nonzero_real = maximum(real.(nonzero_eigs))

    return (
        q6 = pars.q6,
        q5 = pars.q5,
        b6 = 10.0^pars.q6,
        b5 = 10.0^pars.q5,
        T = x[1],
        M = x[2],
        NK = x[3],
        CTL = x[4],
        B = x[5],
        logT = log10(x[1] + 1.0),
        residual = residual,
        smallest_singular_value = S.S[end],
        second_smallest_singular_value = S.S[end-1],
        b20 = b20,
        leading_nonzero_real_eigenvalue = leading_nonzero_real
    )
end

println("\nContinuing the saddle-node in the b5-b6 plane...")

fold_curve = continuation(
    br,
    indfold,
    (@optic _.q5),
    two_parameter_options;
    # Use the default PALC tangent for the fold continuation rather than
    # inheriting the Bordered tangent used for the one-parameter branch.
    alg=PALC(),
    start_with_eigen=true,
    detect_codim2_bifurcation=2,
    update_minaug_every_step=1,
    bdlinsolver=MatrixBLS(),
    bothside=true,
    record_from_solution=record_fold
)

open(joinpath(outdir, "two_parameter_fold_summary.txt"), "w") do io
    show(io, fold_curve)
end

# ----------------------------
# Export fold-curve points
# ----------------------------
fold_df = DataFrame()

for (i, pt) in enumerate(fold_curve.branch)
    push!(fold_df, (
        index = i,
        step = pt.step,
        q6 = pt.q6,
        q5 = pt.q5,
        b6 = pt.b6,
        b5 = pt.b5,
        T = pt.T,
        M = pt.M,
        NK = pt.NK,
        CTL = pt.CTL,
        B = pt.B,
        logT = pt.logT,
        residual = pt.residual,
        smallest_singular_value = pt.smallest_singular_value,
        second_smallest_singular_value = pt.second_smallest_singular_value,
        b20 = pt.b20,
        leading_nonzero_real_eigenvalue =
            pt.leading_nonzero_real_eigenvalue
    ); cols=:union)
end

sort!(fold_df, :b5)
CSV.write(joinpath(outdir, "saddle_node_curve_b5_b6.csv"), fold_df)

# ----------------------------
# Export codimension-two special points
# ----------------------------
special_df = DataFrame()

for (i, sp) in enumerate(fold_curve.specialpoint)
    q5_sp = Float64(sp.param)
    nearest = fold_df[argmin(abs.(fold_df.q5 .- q5_sp)), :]

    push!(special_df, (
        index = i,
        type = String(sp.type),
        status = String(sp.status),
        q5 = q5_sp,
        b5 = 10.0^q5_sp,
        q6_nearest = nearest.q6,
        b6_nearest = nearest.b6,
        T_nearest = nearest.T,
        continuation_step = sp.step,
        precision = sp.precision
    ); cols=:union)
end

CSV.write(joinpath(outdir, "two_parameter_special_points.csv"), special_df)

# ----------------------------
# Optional empirical Step 3 boundary
# ----------------------------
empirical_df = nothing

if empirical_boundary_file !== nothing
    empirical_df = CSV.read(empirical_boundary_file, DataFrame)
    println("Overlaying empirical boundary from: ", empirical_boundary_file)
end

# ----------------------------
# Main figure: formal fold boundary
# ----------------------------
p_curve = plot(
    fold_df.b5,
    fold_df.b6,
    xscale=:log10,
    yscale=:log10,
    marker=:circle,
    ms=4.5,
    lw=2.8,
    color=:darkorange,
    label="formal saddle-node boundary",
    xlabel="b5",
    ylabel="critical b6",
    title="MDSC suppression shifts the saddle-node threshold for tumor control",
    size=(950, 650),
    legend=:topleft,
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm,
    right_margin=8Plots.mm
)

# Baseline fold point.
scatter!(
    p_curve,
    [10.0^Q5_BASE],
    [10.0^sp0.param],
    marker=:star5,
    ms=11,
    color=:black,
    label="baseline saddle-node"
)

# Earlier simulation-defined strict-control boundary.
if empirical_df !== nothing &&
   (:b5 in propertynames(empirical_df)) &&
   (:b6_control_boundary in propertynames(empirical_df))

    empirical_valid = empirical_df[
        .!ismissing.(empirical_df.b6_control_boundary) .&
        .!isnan.(coalesce.(empirical_df.b6_control_boundary, NaN)),
        :
    ]

    if nrow(empirical_valid) > 0
        plot!(
            p_curve,
            empirical_valid.b5,
            empirical_valid.b6_control_boundary,
            marker=:diamond,
            lw=2.0,
            linestyle=:dash,
            color=:steelblue,
            label="simulation T(10000) < 1e3 boundary"
        )
    end
end

# Mark detected cusp / codimension-two points if present.
for sp_type in unique(special_df.type)
    if sp_type in ["cusp", "bt", "zh"]
        sub = special_df[special_df.type .== sp_type, :]
        scatter!(
            p_curve,
            sub.b5,
            sub.b6_nearest,
            marker=:star5,
            ms=10,
            label=uppercase(sp_type)
        )
    end
end

savefig(
    p_curve,
    joinpath(outdir, "saddle_node_threshold_in_b5_b6_plane.png")
)

# ----------------------------
# Supporting figures
# ----------------------------
p_tumor = plot(
    fold_df.b5,
    fold_df.T,
    xscale=:log10,
    yscale=:log10,
    marker=:circle,
    lw=2.6,
    color=:firebrick,
    xlabel="b5",
    ylabel="Tumor equilibrium at saddle-node",
    title="Tumor burden at the saddle-node changes with MDSC suppression",
    legend=false,
    size=(900, 560),
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)
savefig(p_tumor, joinpath(outdir, "tumor_equilibrium_along_saddle_node_curve.png"))

p_b20 = plot(
    fold_df.b5,
    fold_df.b20,
    xscale=:log10,
    marker=:circle,
    lw=2.6,
    color=:purple,
    xlabel="b5",
    ylabel="quadratic fold coefficient b20",
    title="Quadratic fold coefficient along the saddle-node boundary",
    legend=false,
    size=(900, 560),
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)
hline!(p_b20, [0.0], linestyle=:dash, color=:black)
savefig(p_b20, joinpath(outdir, "quadratic_fold_coefficient_along_boundary.png"))

p_residual = plot(
    fold_df.b5,
    max.(fold_df.residual, 1.0e-16),
    xscale=:log10,
    yscale=:log10,
    marker=:circle,
    lw=2.4,
    color=:gray35,
    xlabel="b5",
    ylabel="equilibrium residual",
    title="Fold-curve equilibria satisfy the steady-state equations",
    legend=false,
    size=(900, 540),
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)
hline!(p_residual, [1.0e-8], linestyle=:dash, color=:black)
savefig(p_residual, joinpath(outdir, "saddle_node_curve_residuals.png"))

# ----------------------------
# Concise numerical summary
# ----------------------------
summary = DataFrame(
    quantity = [
        "number_of_fold_curve_points",
        "minimum_b5",
        "maximum_b5",
        "minimum_critical_b6",
        "maximum_critical_b6",
        "number_of_detected_codim2_special_points",
        "number_of_detected_cusps"
    ],
    value = [
        nrow(fold_df),
        minimum(fold_df.b5),
        maximum(fold_df.b5),
        minimum(fold_df.b6),
        maximum(fold_df.b6),
        nrow(special_df),
        sum(special_df.type .== "cusp")
    ]
)
CSV.write(joinpath(outdir, "saddle_node_curve_summary.csv"), summary)

println("\nTwo-parameter saddle-node continuation complete.")
println("Outputs saved to: ", outdir)
println("Main figure:")
println("  saddle_node_threshold_in_b5_b6_plane.png")
println("Supporting figures:")
println("  tumor_equilibrium_along_saddle_node_curve.png")
println("  quadratic_fold_coefficient_along_boundary.png")
println("  saddle_node_curve_residuals.png")
println("Tables:")
println("  saddle_node_curve_b5_b6.csv")
println("  two_parameter_special_points.csv")
