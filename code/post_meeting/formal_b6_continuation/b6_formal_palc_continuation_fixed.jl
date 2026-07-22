using BifurcationKit
using ForwardDiff
using LinearAlgebra
using CSV
using DataFrames
using Plots

# ------------------------------------------------------------
# Formal equilibrium continuation in b6
# ------------------------------------------------------------
# Goal:
#   Trace equilibrium branches through turning points using
#   pseudo-arclength continuation (PALC).
#
# This analysis tests whether the hysteresis-like result is supported by:
#   - coexisting stable equilibrium branches,
#   - an unstable connecting branch,
#   - one or more detected fold points.
#
# Continuation parameter:
#   q = log10(b6)
#
# State coordinates:
#   y_i = log(x_i / scale_i)
#
# This log-state transformation:
#   - keeps all biological populations positive,
#   - improves numerical conditioning,
#   - preserves equilibrium locations,
#   - preserves local stability eigenvalues at equilibria through
#     a similarity transformation.
#
# Figure titles are descriptive only. No "Phase" or "Step" prefixes.
# ------------------------------------------------------------

outdir = "formal_b6_continuation_outputs"
mkpath(outdir)

# ----------------------------
# User settings
# ----------------------------
B5_FIXED = 1.0e-4

# q = log10(b6). Start with the biologically/numerically tested range.
# Widen only after the continuation is working reliably.
Q_MIN = -7.0
Q_MAX = log10(3.0e-5)

# Existing numerical branches provide accurate equilibrium seeds.
HIGH_SEED_B6 = 5.0e-7
LOW_SEED_B6  = 1.5e-6

# Plotting only; does not change continuation data.
TUMOR_PLOT_FLOOR = 1.0e-6
T_CONTROL = 1.0e3
T_DOMINANT = 1.0e5

# Scales used in y = log(x / SCALE).
# They are chosen near the typical magnitudes of each state.
const SCALE = [1.0e7, 5.0e5, 2.0e5, 1.0e3, 1.0e5]

# ----------------------------
# Locate existing branch CSV
# ----------------------------
function first_existing(paths)
    for p in paths
        if isfile(p)
            return p
        end
    end
    error("Could not find b6_equilibrium_branches.csv.\nChecked:\n" *
          join(paths, "\n"))
end

seed_file = first_existing([
    "b6_equilibrium_branches.csv",
    joinpath("b6_equilibrium_jacobian_outputs", "b6_equilibrium_branches.csv"),
    joinpath("results", "post_meeting", "b6_equilibrium_jacobian",
             "b6_equilibrium_branches.csv"),
    joinpath("results", "post_meeting", "b6_equilibrium_jacobian",
             "b6_equilibrium_jacobian_outputs", "b6_equilibrium_branches.csv")
])

seed_rows = CSV.read(seed_file, DataFrame)
println("Using equilibrium seeds from: ", seed_file)

# ----------------------------
# Model
# ----------------------------
function base_p(; b5=B5_FIXED, b6=1.0e-6)
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

    # All continuation states are positive because x = SCALE .* exp.(y).
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

# Convert between scaled log states and biological populations.
x_from_y(y) = SCALE .* exp.(y)
y_from_x(x) = log.(x ./ SCALE)

# IMPORTANT: do not cast the continuation parameter to Float64 here.
# BifurcationKit/ForwardDiff may pass a Dual number when differentiating
# with respect to q. Converting it to Float64 destroys the derivative and
# can make the PALC bordered system singular/non-finite.
@inline function scalar_q(par)
    if par isa Number
        return par
    elseif hasproperty(par, :q)
        return getproperty(par, :q)
    else
        return par[1]
    end
end

# Transformed ODE: dy/dt = (dx/dt) ./ x.
# At an equilibrium, its Jacobian is similar to the original x-space
# Jacobian, so local stability eigenvalues are preserved.
function F(y, par)
    q = scalar_q(par)
    b6 = exp(log(10.0) * q)  # AD-safe equivalent of 10^q

    all(isfinite, y) || error("Non-finite log-state supplied to F: $y")
    isfinite(q) || error("Non-finite continuation parameter q=$q")

    x = x_from_y(y)
    all(isfinite, x) || error("exp(y) overflowed/underflowed in F; y=$y")

    p = base_p(b6=b6)
    out = rhs_x(x, p) ./ x
    all(isfinite, out) || error("Non-finite residual in F at q=$q, x=$x, F=$out")
    return out
end

function J(y, par)
    JJ = ForwardDiff.jacobian(yy -> F(yy, par), y)
    all(isfinite, JJ) || error("Non-finite Jacobian at q=$(scalar_q(par)), y=$y")
    return JJ
end

function record_solution(y, par; k...)
    q = scalar_q(par)
    x = x_from_y(y)
    return (
        b6 = 10.0^q,
        T = x[1],
        M = x[2],
        NK = x[3],
        CTL = x[4],
        B = x[5],
        logT = log10(x[1] + 1.0)
    )
end

# ----------------------------
# Seed selection
# ----------------------------
function nearest_seed(direction::String, target_b6::Float64)
    sub = seed_rows[seed_rows.direction .== direction, :]
    nrow(sub) > 0 || error("No rows found for direction=$(direction).")

    idx = argmin(abs.(log10.(sub.b6) .- log10(target_b6)))
    row = sub[idx, :]

    x = Float64[row.T, row.M, row.NK, row.CTL, row.B]
    q = log10(Float64(row.b6))

    println("Selected $(direction) seed:")
    println("  requested b6 = ", target_b6)
    println("  actual b6    = ", row.b6)
    println("  tumor        = ", row.T)
    println("  residual     = ", row.residual_scaled)

    return y_from_x(x), q
end

y_high, q_high = nearest_seed("up", HIGH_SEED_B6)
y_low,  q_low  = nearest_seed("down", LOW_SEED_B6)

# ----------------------------
# BifurcationKit problems
# ----------------------------
# Use a NamedTuple parameter and an explicit optic. This is the idiomatic
# BifurcationKit setup and keeps q differentiable.
par_high = (q = q_high,)
par_low  = (q = q_low,)
q_lens = (@optic _.q)

prob_high = ODEBifProblem(
    F,
    y_high,
    par_high,
    q_lens;
    J=J,
    record_from_solution=record_solution
)

prob_low = ODEBifProblem(
    F,
    y_low,
    par_low,
    q_lens;
    J=J,
    record_from_solution=record_solution
)

# Preflight diagnostics: stop here with a useful message instead of an
# opaque LAPACK "matrix contains Infs or NaNs" error.
function preflight(name, y0, q0)
    par0 = (q = q0,)
    f0 = F(y0, par0)
    j0 = J(y0, par0)

    hq = 1.0e-6
    dfdq = (F(y0, (q = q0 + hq,)) - F(y0, (q = q0 - hq,))) / (2hq)

    println("\nPreflight: ", name)
    println("  q = ", q0, ", b6 = ", exp(log(10.0) * q0))
    println("  ||F||∞ = ", norm(f0, Inf))
    println("  all finite J = ", all(isfinite, j0))
    println("  ||dF/dq||∞ = ", norm(dfdq, Inf))
    println("  cond(J) = ", cond(j0))

    all(isfinite, f0) || error("$name seed has non-finite residual")
    all(isfinite, j0) || error("$name seed has non-finite Jacobian")
    all(isfinite, dfdq) || error("$name seed has non-finite dF/dq")
    norm(dfdq, Inf) > 0 || error("$name has zero dF/dq; parameter dependence was lost")
end

preflight("high-tumor seed", y_high, q_high)
preflight("low-tumor seed", y_low, q_low)

newton_options = NewtonPar(
    tol=1.0e-10,
    max_iterations=40,
    linesearch=true,
    verbose=true
)

continuation_options = ContinuationPar(
    ds=1.0e-3,
    dsmin=1.0e-7,
    dsmax=1.0e-2,
    p_min=Q_MIN,
    p_max=Q_MAX,
    max_steps=3000,

    newton_options=newton_options,

    # Five-dimensional system: compute all five eigenvalues.
    nev=5,
    save_eig_every_step=1,
    save_eigenvectors=false,

    # Detect folds and stability changes precisely.
    detect_fold=true,
    detect_bifurcation=3,
    n_inversion=6,
    max_bisection_steps=30,
    tol_bisection_eigenvalue=1.0e-8,
    tol_stability=1.0e-8
)

println("\nRunning PALC from the high-tumor equilibrium seed...")
branch_high = continuation(
    prob_high,
    PALC(tangent=Secant(), θ=0.1),
    continuation_options;
    bothside=true,
    verbosity=3
)

println("\nRunning PALC from the low-tumor equilibrium seed...")
branch_low = continuation(
    prob_low,
    PALC(tangent=Secant(), θ=0.1),
    continuation_options;
    bothside=true,
    verbosity=3
)

# Save BifurcationKit text summaries.
open(joinpath(outdir, "high_seed_continuation_summary.txt"), "w") do io
    show(io, branch_high)
end

open(joinpath(outdir, "low_seed_continuation_summary.txt"), "w") do io
    show(io, branch_low)
end

# ----------------------------
# Export branch points
# ----------------------------
function branch_dataframe(br, seed_label)
    df = DataFrame()

    for (i, pt) in enumerate(br.branch)
        q = scalar_q(pt.param)
        b6 = 10.0^q

        x = Float64[pt.T, pt.M, pt.NK, pt.CTL, pt.B]
        y = y_from_x(x)

        residual = norm(F(y, [q]), Inf)
        eigs = eigvals(J(y, [q]))
        max_real_eig = maximum(real.(eigs))
        stable = max_real_eig < 0.0

        push!(df, (
            seed = seed_label,
            branch_index = i,
            step = pt.step,
            q = q,
            b6 = b6,
            T = x[1],
            M = x[2],
            NK = x[3],
            CTL = x[4],
            B = x[5],
            logT = log10(x[1] + 1.0),
            residual = residual,
            max_real_eig = max_real_eig,
            stable = stable,
            n_unstable = pt.n_unstable,
            n_imag = pt.n_imag
        ); cols=:union)
    end

    return df
end

high_df = branch_dataframe(branch_high, "high_seed")
low_df  = branch_dataframe(branch_low, "low_seed")
all_df  = vcat(high_df, low_df; cols=:union)

CSV.write(joinpath(outdir, "formal_continuation_high_seed.csv"), high_df)
CSV.write(joinpath(outdir, "formal_continuation_low_seed.csv"), low_df)
CSV.write(joinpath(outdir, "formal_continuation_all_points.csv"), all_df)

# ----------------------------
# Export special points
# ----------------------------
function getfield_or(x, name::Symbol, default)
    return hasproperty(x, name) ? getproperty(x, name) : default
end

function specialpoint_dataframe(br, seed_label)
    df = DataFrame()

    for sp in br.specialpoint
        q = scalar_q(sp.param)
        psol = getfield_or(sp, :printsol, nothing)
        T = (psol !== nothing && hasproperty(psol, :T)) ? psol.T : NaN

        push!(df, (
            seed = seed_label,
            type = String(sp.type),
            status = String(sp.status),
            q = q,
            b6 = 10.0^q,
            T = T,
            branch_index = sp.idx,
            continuation_step = sp.step,
            precision = sp.precision,
            interval_left = sp.interval[1],
            interval_right = sp.interval[2],
            delta_real = sp.δ[1],
            delta_imag = sp.δ[2]
        ); cols=:union)
    end

    return df
end

sp_high = specialpoint_dataframe(branch_high, "high_seed")
sp_low  = specialpoint_dataframe(branch_low, "low_seed")
sp_all  = vcat(sp_high, sp_low; cols=:union)
CSV.write(joinpath(outdir, "formal_continuation_special_points.csv"), sp_all)

# ----------------------------
# Plot helpers
# ----------------------------
plot_T(T) = max.(T, TUMOR_PLOT_FLOOR)

function add_branch!(p, df; stable_color, unstable_color,
                     stable_label, unstable_label)
    # A faint line shows branch ordering through any folds.
    plot!(p, df.b6, plot_T(df.T),
          color=:gray, alpha=0.35, lw=1.2, label=false)

    stable_df = df[df.stable .== true, :]
    unstable_df = df[df.stable .== false, :]

    if nrow(stable_df) > 0
        scatter!(p, stable_df.b6, plot_T(stable_df.T),
                 marker=:circle, ms=4.5,
                 color=stable_color, label=stable_label)
    end

    if nrow(unstable_df) > 0
        scatter!(p, unstable_df.b6, plot_T(unstable_df.T),
                 marker=:circle, ms=5,
                 markerstrokecolor=unstable_color,
                 markercolor=:white,
                 label=unstable_label)
    end
end

# ----------------------------
# Figure 1: equilibrium branches
# ----------------------------
p_branch = plot(
    xscale=:log10,
    yscale=:log10,
    xlabel="b6",
    ylabel="Tumor equilibrium",
    title="Equilibrium branches across b6",
    size=(950, 650),
    legend=:bottomleft,
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm,
    right_margin=8Plots.mm
)

add_branch!(
    p_branch, high_df;
    stable_color=:steelblue,
    unstable_color=:navy,
    stable_label="high-seed stable",
    unstable_label="high-seed unstable"
)

add_branch!(
    p_branch, low_df;
    stable_color=:darkorange,
    unstable_color=:red,
    stable_label="low-seed stable",
    unstable_label="low-seed unstable"
)

hline!(p_branch, [T_CONTROL],
       linestyle=:dash, color=:black, label="T = 1e3")
hline!(p_branch, [T_DOMINANT],
       linestyle=:dot, color=:gray, label="T = 1e5")

# Mark detected folds.
folds = sp_all[sp_all.type .== "fold", :]
if nrow(folds) > 0
    scatter!(p_branch, folds.b6, plot_T(folds.T),
             marker=:star5, ms=10,
             color=:red, label="detected fold")
end

savefig(p_branch, joinpath(outdir, "equilibrium_branches_across_b6.png"))

# ----------------------------
# Figure 2: local stability
# ----------------------------
p_stability = plot(
    high_df.b6, high_df.max_real_eig,
    xscale=:log10,
    marker=:circle, ms=3.5, lw=2,
    color=:steelblue,
    label="high-seed branch",
    xlabel="b6",
    ylabel="maximum real eigenvalue",
    title="Local stability along equilibrium branches",
    size=(950, 570),
    legend=:bottomright,
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)

plot!(p_stability, low_df.b6, low_df.max_real_eig,
      marker=:diamond, ms=3.5, lw=2,
      color=:darkorange,
      label="low-seed branch")

hline!(p_stability, [0.0],
       linestyle=:dash, color=:black,
       label="stability boundary")

savefig(p_stability, joinpath(outdir, "local_stability_along_branches.png"))

# ----------------------------
# Figure 3: residual quality
# ----------------------------
safe_residual(r) = max.(r, 1.0e-16)

p_residual = plot(
    high_df.b6, safe_residual(high_df.residual),
    xscale=:log10,
    yscale=:log10,
    marker=:circle, ms=3.5, lw=2,
    color=:steelblue,
    label="high-seed branch",
    xlabel="b6",
    ylabel="||F(y,q)||∞",
    title="Equilibrium residuals along continuation branches",
    size=(950, 570),
    legend=:bottomleft,
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)

plot!(p_residual, low_df.b6, safe_residual(low_df.residual),
      marker=:diamond, ms=3.5, lw=2,
      color=:darkorange,
      label="low-seed branch")

hline!(p_residual, [1.0e-8],
       linestyle=:dash, color=:black,
       label="1e-8")

savefig(p_residual, joinpath(outdir, "equilibrium_residuals_along_branches.png"))

# ----------------------------
# Summary text
# ----------------------------
open(joinpath(outdir, "formal_continuation_readme.txt"), "w") do io
    println(io, "Formal b6 equilibrium continuation")
    println(io, "=================================")
    println(io)
    println(io, "Continuation variable: q = log10(b6)")
    println(io, "Parameter range: b6 = 1e$(Q_MIN) to 1e$(Q_MAX)")
    println(io, "High-tumor seed target b6: ", HIGH_SEED_B6)
    println(io, "Low-tumor seed target b6: ", LOW_SEED_B6)
    println(io)
    println(io, "High-seed points: ", nrow(high_df))
    println(io, "Low-seed points: ", nrow(low_df))
    println(io, "Detected special points: ", nrow(sp_all))
    println(io, "Detected folds: ", sum(sp_all.type .== "fold"))
    println(io)
    println(io, "Interpretation rule:")
    println(io, "- stable points have max real eigenvalue < 0")
    println(io, "- unstable points have max real eigenvalue > 0")
    println(io, "- fold points are candidate saddle-node transitions")
    println(io, "- a stable–unstable–stable branch structure would support bistability")
end

println("\nFormal continuation complete.")
println("Outputs saved to: ", outdir)
println("Main figure: equilibrium_branches_across_b6.png")
println("Other figures:")
println("  local_stability_along_branches.png")
println("  equilibrium_residuals_along_branches.png")
println("Tables:")
println("  formal_continuation_all_points.csv")
println("  formal_continuation_special_points.csv")
