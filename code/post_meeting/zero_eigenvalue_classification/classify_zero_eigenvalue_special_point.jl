using LinearAlgebra
using ForwardDiff
using CSV
using DataFrames
using Plots
using Printf

# ------------------------------------------------------------
# Local classification of the zero-eigenvalue equilibrium
# ------------------------------------------------------------
# Purpose:
#   Refine the special point detected on the low-tumor equilibrium branch and
#   test whether it satisfies the generic saddle-node (fold) conditions.
#
# Professional terminology:
#   - equilibrium / fixed point: F(x, b6) = 0
#   - nonhyperbolic equilibrium: Jacobian has an eigenvalue at zero
#   - saddle-node / fold bifurcation: a generic codimension-one
#     zero-eigenvalue equilibrium where the equilibrium branch turns
#
# The script tests:
#   1. equilibrium residual is small
#   2. the Jacobian has one isolated zero singular value/eigenvalue
#   3. a01 = w' * F_q is nonzero (parameter transversality)
#   4. b20 = w' * D²F[v,v] is nonzero (quadratic nondegeneracy)
#
# For a simple zero eigenvalue:
#   a01 ≠ 0 and b20 ≠ 0  => generic saddle-node / fold
#
# Figure titles are descriptive only; no workflow prefixes.
# ------------------------------------------------------------

outdir = "zero_eigenvalue_classification_outputs"
mkpath(outdir)

# ----------------------------
# User settings
# ----------------------------
B5_FIXED = 1.0e-4

# Numerical tolerances used only for classification.
SIMPLE_ZERO_RATIO_TOL = 1.0e-5
COEFFICIENT_TOL = 1.0e-6
NEWTON_TOL = 1.0e-11
MAX_NEWTON_ITER = 40

# Scaling from the formal continuation analysis.
const SCALE = [1.0e7, 5.0e5, 2.0e5, 1.0e3, 1.0e5]

# ----------------------------
# Locate input files
# ----------------------------
function first_existing(paths)
    for p in paths
        if isfile(p)
            return p
        end
    end
    error("Could not find any expected file:\n" * join(paths, "\n"))
end

low_file = first_existing([
    "formal_continuation_low_seed.csv",
    joinpath("formal_b6_continuation_outputs", "formal_continuation_low_seed.csv"),
    joinpath("results", "post_meeting", "formal_b6_continuation",
             "formal_continuation_low_seed.csv")
])

special_file = first_existing([
    "formal_continuation_special_points.csv",
    joinpath("formal_b6_continuation_outputs", "formal_continuation_special_points.csv"),
    joinpath("results", "post_meeting", "formal_b6_continuation",
             "formal_continuation_special_points.csv")
])

low = CSV.read(low_file, DataFrame)
special = CSV.read(special_file, DataFrame)

println("Using low branch: ", low_file)
println("Using special points: ", special_file)

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

# q = log10(b6)
function F(y, q)
    x = x_from_y(y)
    p = base_p(b6=10.0^q)
    return rhs_x(x, p) ./ x
end

function Jy(y, q)
    return ForwardDiff.jacobian(yy -> F(yy, q), y)
end

# ----------------------------
# Find initial special-point estimate
# ----------------------------
candidate_rows = special[special.type .== "bp", :]

if nrow(candidate_rows) > 0
    q_guess = candidate_rows.q[1]
else
    # Fallback: use branch point with eigenvalue closest to zero.
    idx0 = argmin(abs.(low.max_real_eig))
    q_guess = low.q[idx0]
end

idx_near = argmin(abs.(low.q .- q_guess))
row0 = low[idx_near, :]

x0 = Float64[row0.T, row0.M, row0.NK, row0.CTL, row0.B]
y0 = y_from_x(x0)
q0 = Float64(row0.q)

J0 = Jy(y0, q0)
sv0 = svd(J0)
v0 = sv0.V[:, end]
v0 ./= norm(v0)

println("\nInitial zero-eigenvalue candidate:")
println("  q = ", q0)
println("  b6 = ", 10.0^q0)
println("  tumor = ", x0[1])
println("  max real eigenvalue = ", row0.max_real_eig)
println("  continuation residual = ", row0.residual)

# ----------------------------
# Augmented fold equations
# ----------------------------
# Unknown z = [y(1:5), q, v(1:5)]
#
# Equations:
#   F(y,q) = 0
#   J(y,q) v = 0
#   ||v||² - 1 = 0
#
# This directly refines a zero-eigenvalue equilibrium.
function fold_residual(z)
    y = z[1:5]
    q = z[6]
    v = z[7:11]

    J = Jy(y, q)

    return vcat(
        F(y, q),
        J * v,
        dot(v, v) - 1.0
    )
end

function fd_jacobian(fun, z; relstep=1.0e-6)
    m = length(fun(z))
    n = length(z)
    A = zeros(m, n)

    for j in 1:n
        h = relstep * max(abs(z[j]), 1.0)
        zp = copy(z)
        zm = copy(z)
        zp[j] += h
        zm[j] -= h
        A[:, j] = (fun(zp) - fun(zm)) / (2h)
    end

    return A
end

function refine_fold(z0; tol=NEWTON_TOL, maxiter=MAX_NEWTON_ITER)
    z = copy(z0)
    best_z = copy(z)
    best_norm = norm(fold_residual(z), Inf)
    converged = false

    for iter in 1:maxiter
        r = fold_residual(z)
        rn = norm(r, Inf)

        @printf("Newton iteration %2d: residual = %.4e\n", iter, rn)

        if rn < best_norm
            best_norm = rn
            best_z = copy(z)
        end

        if rn < tol
            converged = true
            best_z = copy(z)
            best_norm = rn
            break
        end

        A = fd_jacobian(fold_residual, z)

        delta = try
            -(A + 1.0e-12I) \ r
        catch
            -(A + 1.0e-8I) \ r
        end

        # Prevent very large Newton jumps.
        maxstep = maximum(abs.(delta))
        if maxstep > 1.0
            delta .*= 1.0 / maxstep
        end

        # Backtracking line search.
        alpha = 1.0
        accepted = false

        for _ in 1:16
            z_try = z + alpha * delta
            rn_try = norm(fold_residual(z_try), Inf)

            if isfinite(rn_try) && rn_try < rn
                z = z_try
                accepted = true
                break
            end

            alpha *= 0.5
        end

        if !accepted
            println("Line search stalled.")
            break
        end
    end

    return best_z, converged, best_norm
end

z0 = vcat(y0, q0, v0)
zstar, converged, augmented_residual = refine_fold(z0)

ystar = zstar[1:5]
qstar = zstar[6]
xstar = x_from_y(ystar)
b6star = 10.0^qstar

# ----------------------------
# Local nullspace and normal-form coefficients
# ----------------------------
Jstar = Jy(ystar, qstar)
sv = svd(Jstar)

# Right and left nullvectors from the smallest singular value.
v = sv.V[:, end]
w = sv.U[:, end]

# Normalize so w'v = 1.
if dot(w, v) < 0
    w .*= -1.0
end
w ./= dot(w, v)
v ./= norm(v)

eigs = eigvals(Jstar)
idx_eig = sortperm(abs.(eigs))

# Parameter derivative F_q.
hq = 1.0e-6
Fq = (F(ystar, qstar + hq) - F(ystar, qstar - hq)) / (2hq)

# Quadratic directional derivative D²F[v,v].
hv = 1.0e-4
Hvv = (
    F(ystar + hv*v, qstar) -
    2.0F(ystar, qstar) +
    F(ystar - hv*v, qstar)
) / hv^2

# Cubic directional derivative D³F[v,v,v], retained as a diagnostic.
D3vvv = (
    F(ystar + 2hv*v, qstar) -
    2.0F(ystar + hv*v, qstar) +
    2.0F(ystar - hv*v, qstar) -
    F(ystar - 2hv*v, qstar)
) / (2hv^3)

# Standard fold nondegeneracy diagnostics in these scaled coordinates.
a01 = dot(w, Fq)
b20 = dot(w, Hvv)
b30 = dot(w, D3vvv)

simple_zero_ratio = sv.S[end] / sv.S[end-1]
simple_zero = simple_zero_ratio < SIMPLE_ZERO_RATIO_TOL

classification = if simple_zero &&
                    abs(a01) > COEFFICIENT_TOL &&
                    abs(b20) > COEFFICIENT_TOL
    "generic saddle-node (fold) bifurcation"
elseif simple_zero && abs(a01) <= COEFFICIENT_TOL
    "branch-point candidate; compute full Lyapunov-Schmidt normal form"
else
    "degenerate or unresolved zero-eigenvalue point"
end

# ----------------------------
# Save numerical results
# ----------------------------
summary = DataFrame(
    quantity = [
        "classification",
        "newton_converged",
        "augmented_residual",
        "q_star",
        "b6_star",
        "Tumor_star",
        "MDSC_star",
        "NK_star",
        "CTL_star",
        "B_star",
        "smallest_singular_value",
        "second_smallest_singular_value",
        "simple_zero_ratio",
        "a01_parameter_transversality",
        "b20_quadratic_nondegeneracy",
        "b30_cubic_diagnostic"
    ],
    value = Any[
        classification,
        converged,
        augmented_residual,
        qstar,
        b6star,
        xstar[1],
        xstar[2],
        xstar[3],
        xstar[4],
        xstar[5],
        sv.S[end],
        sv.S[end-1],
        simple_zero_ratio,
        a01,
        b20,
        b30
    ]
)

CSV.write(joinpath(outdir, "zero_eigenvalue_classification.csv"), summary)

eig_df = DataFrame(
    real_part = real.(eigs),
    imaginary_part = imag.(eigs),
    magnitude = abs.(eigs)
)
sort!(eig_df, :magnitude)
CSV.write(joinpath(outdir, "zero_eigenvalue_eigenvalues.csv"), eig_df)

open(joinpath(outdir, "zero_eigenvalue_classification.txt"), "w") do io
    println(io, "Zero-eigenvalue equilibrium classification")
    println(io, "==========================================")
    println(io)
    println(io, "Classification: ", classification)
    println(io, "Converged: ", converged)
    println(io, "Augmented residual: ", augmented_residual)
    println(io)
    println(io, "Critical equilibrium:")
    println(io, "  q* = ", qstar)
    println(io, "  b6* = ", b6star)
    println(io, "  Tumor* = ", xstar[1])
    println(io, "  MDSC* = ", xstar[2])
    println(io, "  NK* = ", xstar[3])
    println(io, "  CTL* = ", xstar[4])
    println(io, "  B* = ", xstar[5])
    println(io)
    println(io, "Simple-zero diagnostics:")
    println(io, "  smallest singular value = ", sv.S[end])
    println(io, "  second-smallest singular value = ", sv.S[end-1])
    println(io, "  ratio = ", simple_zero_ratio)
    println(io)
    println(io, "Nondegeneracy coefficients:")
    println(io, "  a01 = w'F_q = ", a01)
    println(io, "  b20 = w'D²F[v,v] = ", b20)
    println(io, "  b30 = w'D³F[v,v,v] = ", b30)
    println(io)
    println(io, "Interpretation:")
    println(io, "  A simple zero eigenvalue plus nonzero a01 and b20")
    println(io, "  satisfies the generic local conditions for a saddle-node/fold.")
end

# ----------------------------
# Local figures
# ----------------------------
# Select a local continuation window around the refined point.
local_mask = abs.(low.q .- qstar) .< 0.09
local_branch = low[local_mask, :]
sort!(local_branch, :branch_index)

stable_local = local_branch[local_branch.max_real_eig .< 0.0, :]
unstable_local = local_branch[local_branch.max_real_eig .> 0.0, :]

# 1. Branch geometry
p_branch = plot(
    xlabel="b6",
    ylabel="Tumor equilibrium",
    xscale=:log10,
    yscale=:log10,
    title=(occursin("saddle-node", classification) ?
           "Stable and unstable tumor equilibria meet at a saddle-node" :
           "Stable and unstable tumor equilibria meet near a zero-eigenvalue point"),
    size=(900, 600),
    legend=:bottomleft,
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)

if nrow(stable_local) > 0
    plot!(p_branch, stable_local.b6, stable_local.T,
          marker=:circle, lw=2.8, color=:darkorange,
          label="stable branch")
end

if nrow(unstable_local) > 0
    plot!(p_branch, unstable_local.b6, unstable_local.T,
          marker=:circle, markercolor=:white,
          markerstrokecolor=:red, lw=2.8, color=:red,
          label="unstable branch")
end

scatter!(p_branch, [b6star], [xstar[1]],
         marker=:star5, ms=11, color=:black,
         label="zero-eigenvalue equilibrium")

savefig(p_branch, joinpath(outdir, "stable_and_unstable_equilibria_meet.png"))

# 2. Leading eigenvalue
p_eig = plot(
    local_branch.b6, local_branch.max_real_eig,
    xscale=:log10,
    marker=:circle, lw=2.8,
    color=:purple,
    xlabel="b6",
    ylabel="maximum real eigenvalue",
    title="A simple real eigenvalue crosses zero at the critical b6",
    legend=false,
    size=(900, 550),
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)

hline!(p_eig, [0.0], linestyle=:dash, color=:black)
vline!(p_eig, [b6star], linestyle=:dashdot, color=:red)

savefig(p_eig, joinpath(outdir, "real_eigenvalue_crosses_zero.png"))

# 3. Singular-value spectrum
p_sv = scatter(
    1:length(sv.S), sort(sv.S),
    yscale=:log10,
    marker=:circle, ms=7,
    color=:steelblue,
    xlabel="ordered singular-value index",
    ylabel="singular value of Jacobian",
    title="The Jacobian has one isolated zero mode",
    legend=false,
    size=(800, 520),
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)

savefig(p_sv, joinpath(outdir, "jacobian_has_one_zero_mode.png"))

println("\nClassification complete.")
println("Classification: ", classification)
println("b6* = ", b6star)
println("Tumor* = ", xstar[1])
println("a01 = ", a01)
println("b20 = ", b20)
println("Outputs saved to: ", outdir)
