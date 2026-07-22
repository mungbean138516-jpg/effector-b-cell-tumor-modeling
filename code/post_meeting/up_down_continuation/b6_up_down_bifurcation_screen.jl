using DifferentialEquations
using Plots
using DataFrames
using CSV
using LinearAlgebra

# ------------------------------------------------------------
# Step 4: b6 up-sweep / down-sweep bifurcation screen
# ------------------------------------------------------------
# Goal:
#   Test whether the long-time b6 threshold behaves like a simple sharp switch
#   or shows hysteresis / possible bistability.
#
# Logic:
#   Up-sweep:
#       start at low b6, simulate to long time, then use the final state as
#       the initial condition for the next larger b6.
#   Down-sweep:
#       start at high b6, simulate to long time, then use the final state as
#       the initial condition for the next smaller b6.
#
# Interpretation:
#   - If up and down curves overlap: likely a single stable long-time branch
#     / sharp switch / no clear hysteresis in this tested range.
#   - If up and down curves differ: possible hysteresis / bistability, which
#     motivates a more formal bifurcation analysis.
#
# This is not a formal proof of bifurcation. It is a numerical screening test.
# ------------------------------------------------------------

outdir = "b6_bifurcation_outputs"
mkpath(outdir)

t_control = 1.0e3

default_tmax = 10000.0

function logspace10(lo, hi, n)
    return 10 .^ range(log10(lo), log10(hi), length=n)
end

safe(v; eps=1e-12) = max.(collect(v), eps)

function base_p(; b5=1.0e-4, b6=1.0e-7)
    return (
        # original CIR-style parameters
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

        # effector B extension
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

function baseline_state(p; tumor_cells=2.0)
    xM0 = p.a2 / p.z2
    xN0 = p.z2 * p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC0 = 0.0
    xB0 = p.a8 / (p.z5 + p.b5 * xM0)
    return [tumor_cells, xM0, xN0, xC0, xB0]
end

function model!(du, u, p, t)
    xT, xM, xN, xC, xB = max.(u, 0.0)

    # direct B -> CTL term, fixed from previous analyses
    b_to_ctl = p.a11 * xB / (p.g5 + xB)

    growth = p.a1 * xT * log(max(p.h / max(xT, 1e-12), 1.0))
    du[1] = growth - p.b1*xT*xN - p.b2*xT*xC - p.b6*xT*xB - p.z1*xT
    du[2] = p.a2 + p.a3*xT/(p.g1 + xT) - p.z2*xM

    nk_recruit = p.a5 * xT^2/(p.g2 + xT^2) * (1.0 + p.a10*xB/(p.g5 + xB))
    du[3] = p.a4 + nk_recruit - p.b3*xM*xN - p.z3*xN

    du[4] = p.a6*xT*xN + p.a7*xT^2/(p.g3 + xT^2) + b_to_ctl - p.b4*xM*xC - p.z4*xC
    du[5] = p.a8 + p.a9*xT^2/(p.g4 + xT^2) - p.b5*xM*xB - p.z5*xB
end

function fvec(u, p)
    du = similar(u)
    model!(du, u, p, 0.0)
    return du
end

function residual_norm(u, p)
    return norm(fvec(u, p))
end

# Finite-difference Jacobian at a numerical steady state.
# This is a basic stability screen; it should not replace formal continuation.
function jacobian_fd(u, p; relstep=1e-6)
    n = length(u)
    J = zeros(n, n)
    f0 = fvec(u, p)
    for j in 1:n
        h = relstep * max(abs(u[j]), 1.0)
        up = copy(u)
        up[j] += h
        fp = fvec(up, p)
        J[:, j] = (fp - f0) / h
    end
    return J
end

function max_real_eig(u, p)
    J = jacobian_fd(u, p)
    vals = eigvals(J)
    return maximum(real.(vals))
end

function bad_state(u, t, integrator)
    any(x -> !isfinite(x) || x < -1e-6, u)
end
cb_bad = DiscreteCallback(bad_state, terminate!)

function solve_to_final(u0, b6; b5=1.0e-4, tmax=default_tmax)
    p = base_p(b5=b5, b6=b6)
    prob = ODEProblem(model!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5();
                saveat=[tmax],
                reltol=1e-8,
                abstol=1e-10,
                callback=cb_bad,
                isoutofdomain=(u,p,t)->any(x -> x < -1e-6, u))
    u_final = max.(Array(sol.u[end]), 0.0)
    return u_final, p
end

function classify_state(T)
    if T < 1e3
        return "controlled"
    elseif T < 1e5
        return "intermediate"
    else
        return "tumor_dominant"
    end
end

function continuation_sweep(b6_values; direction="up", b5=1.0e-4, tmax=default_tmax)
    p_start = base_p(b5=b5, b6=b6_values[1])
    u_current = baseline_state(p_start)

    rows = DataFrame(
        direction=String[], b5=Float64[], b6=Float64[],
        T=Float64[], MDSC=Float64[], NK=Float64[], CTL=Float64[], B=Float64[],
        logT=Float64[], state=String[], residual=Float64[], max_real_eig=Float64[]
    )

    for (k, b6) in enumerate(b6_values)
        println(direction, " sweep: b6=", b6, " (", k, "/", length(b6_values), ")")
        u_final, p = solve_to_final(u_current, b6; b5=b5, tmax=tmax)
        res = residual_norm(u_final, p)
        maxeig = max_real_eig(u_final, p)
        push!(rows, (
            direction, b5, b6,
            u_final[1], u_final[2], u_final[3], u_final[4], u_final[5],
            log10(u_final[1] + 1.0), classify_state(u_final[1]), res, maxeig
        ))
        # Continuation step: final state becomes the next initial condition.
        u_current = u_final
    end
    return rows
end

# ------------------------------------------------------------
# Run up-sweep and down-sweep
# ------------------------------------------------------------
# Focus near the long-time threshold found in Step 1.
b6_grid = logspace10(2.5e-7, 2.0e-6, 80)

println("Step 4: up-sweep / down-sweep test")
println("  b5 fixed at 1e-4")
println("  tmax per step = ", default_tmax)
println("  b6 grid points = ", length(b6_grid))

up = continuation_sweep(b6_grid; direction="up", b5=1e-4, tmax=default_tmax)
down = continuation_sweep(reverse(b6_grid); direction="down", b5=1e-4, tmax=default_tmax)

allrows = vcat(up, down)
CSV.write(joinpath(outdir, "b6_up_down_continuation.csv"), allrows)
println("Saved: ", joinpath(outdir, "b6_up_down_continuation.csv"))

# Comparison table at matching b6 values.
# Down sweep is stored descending, so sort before merge by b6.
up_s = sort(up, :b6)
down_s = sort(down, :b6)
comp = DataFrame(
    b6 = up_s.b6,
    T_up = up_s.T,
    T_down = down_s.T,
    logT_up = up_s.logT,
    logT_down = down_s.logT,
    abs_logT_difference = abs.(up_s.logT .- down_s.logT),
    state_up = up_s.state,
    state_down = down_s.state,
    maxeig_up = up_s.max_real_eig,
    maxeig_down = down_s.max_real_eig,
    residual_up = up_s.residual,
    residual_down = down_s.residual,
)
CSV.write(joinpath(outdir, "b6_up_down_comparison.csv"), comp)
println("Saved: ", joinpath(outdir, "b6_up_down_comparison.csv"))

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------
# Main hysteresis plot: T vs b6 for up vs down.
p1 = plot(
    up_s.b6, max.(up_s.T, 1e-8),
    xscale=:log10, yscale=:log10,
    marker=:circle, lw=2.5,
    label="up-sweep",
    xlabel="b6",
    ylabel="T after continuation step",
    title="Up-sweep vs down-sweep: tumor branch",
    size=(900, 600)
)
plot!(p1, down_s.b6, max.(down_s.T, 1e-8),
      marker=:diamond, lw=2.5, label="down-sweep")
hline!(p1, [1e3], linestyle=:dash, color=:black, label="T=1e3")
hline!(p1, [1e5], linestyle=:dot, color=:gray, label="T=1e5")
savefig(p1, joinpath(outdir, "b6_up_down_tumor_branch.png"))

# Difference plot: if nearly zero, no clear hysteresis.
p2 = plot(
    comp.b6, comp.abs_logT_difference,
    xscale=:log10,
    marker=:circle, lw=2.5,
    xlabel="b6",
    ylabel="|logT_up - logT_down|",
    title="Up/down difference: hysteresis screen",
    legend=false,
    size=(850, 550)
)
savefig(p2, joinpath(outdir, "b6_up_down_difference.png"))

# Max real eigenvalue plot: negative suggests local stability of the numerical branch.
p3 = plot(
    up_s.b6, up_s.max_real_eig,
    xscale=:log10,
    marker=:circle, lw=2.5,
    label="up-sweep",
    xlabel="b6",
    ylabel="max real eigenvalue",
    title="Local stability screen near numerical long-time states",
    size=(900, 600)
)
plot!(p3, down_s.b6, down_s.max_real_eig,
      marker=:diamond, lw=2.5, label="down-sweep")
hline!(p3, [0.0], linestyle=:dash, color=:black, label="0")
savefig(p3, joinpath(outdir, "b6_up_down_stability_screen.png"))

println("Step 4 up/down screen complete.")
println("Main files:")
println("  ", joinpath(outdir, "b6_up_down_tumor_branch.png"))
println("  ", joinpath(outdir, "b6_up_down_difference.png"))
println("  ", joinpath(outdir, "b6_up_down_stability_screen.png"))
println("Interpretation guide:")
println("  If up/down tumor curves overlap -> no clear hysteresis in this range.")
println("  If up/down tumor curves differ -> possible bistability / hysteresis.")
