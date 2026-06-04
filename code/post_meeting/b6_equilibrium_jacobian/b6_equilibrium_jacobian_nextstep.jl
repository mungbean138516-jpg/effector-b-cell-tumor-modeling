using DifferentialEquations
using LinearAlgebra
using Statistics
using DataFrames
using CSV
using Plots

# Step 5: steady-state and Jacobian validation for b6.
#
# The up/down b6 sweeps suggested hysteresis, so this script treats the
# long-time states as equilibrium candidates, refines them, and checks local
# stability with finite-difference Jacobians.

outdir = "b6_equilibrium_jacobian_outputs"
mkpath(outdir)

# Utility helpers
function logspace10(lo, hi, n)
    return 10 .^ range(log10(lo), log10(hi), length=n)
end

safe(v; eps=1e-12) = max.(collect(v), eps)

function state_label(T)
    if T < 1e3
        return "controlled"
    elseif T < 1e5
        return "intermediate"
    else
        return "tumor_dominant"
    end
end

# Model parameters
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

# ODE right-hand side
function rhs_vec(x, p)
    # keep states non-negative for function evaluation
    xT, xM, xN, xC, xB = max.(x, 0.0)
    du = zeros(5)

    b_to_ctl = p.a11 * xB / (p.g5 + xB)
    growth = p.a1 * xT * log(max(p.h / max(xT, 1e-12), 1.0))

    du[1] = growth - p.b1*xT*xN - p.b2*xT*xC - p.b6*xT*xB - p.z1*xT
    du[2] = p.a2 + p.a3*xT/(p.g1 + xT) - p.z2*xM

    nk_recruit = p.a5 * xT^2/(p.g2 + xT^2) * (1.0 + p.a10*xB/(p.g5 + xB))
    du[3] = p.a4 + nk_recruit - p.b3*xM*xN - p.z3*xN

    du[4] = p.a6*xT*xN + p.a7*xT^2/(p.g3 + xT^2) + b_to_ctl - p.b4*xM*xC - p.z4*xC
    du[5] = p.a8 + p.a9*xT^2/(p.g4 + xT^2) - p.b5*xM*xB - p.z5*xB
    return du
end

function model!(du, u, p, t)
    du[:] = rhs_vec(u, p)
end

function bad_state(u, t, integrator)
    any(x -> !isfinite(x) || x < -1e-6, u)
end
cb_bad = DiscreteCallback(bad_state, terminate!)

function solve_long(b6, u0; b5=1e-4, tmax=10000.0, saveat=tmax)
    p = base_p(b5=b5, b6=b6)
    prob = ODEProblem(model!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5();
                saveat=[tmax],
                reltol=1e-7,
                abstol=1e-9,
                callback=cb_bad,
                isoutofdomain=(u,p,t)->any(x -> x < -1e-6, u))
    return sol.u[end], p
end

# Finite-difference Jacobian in x-space
function jacobian_fd(x, p; rel_step=1e-5)
    n = length(x)
    J = zeros(n,n)
    f0 = rhs_vec(x, p)
    for j in 1:n
        step = rel_step * max(abs(x[j]), 1.0)
        xp = copy(x); xm = copy(x)
        xp[j] += step
        xm[j] = max(xm[j] - step, 1e-12)
        fp = rhs_vec(xp, p)
        fm = rhs_vec(xm, p)
        J[:,j] = (fp - fm) / (xp[j] - xm[j])
    end
    return J
end

function max_real_eig(x, p)
    vals = eigvals(jacobian_fd(x, p))
    return maximum(real.(vals))
end

# Newton refinement in log-state variables keeps candidate states positive
# and improves conditioning across cell populations.
function scaled_residual_y(y, p)
    x = exp.(y)
    f = rhs_vec(x, p)
    # scale each equation by the typical size of the state + 1
    return f ./ max.(abs.(x), 1.0)
end

function jacobian_y_fd(y, p; rel_step=1e-5)
    n = length(y)
    J = zeros(n,n)
    for j in 1:n
        step = rel_step
        yp = copy(y); ym = copy(y)
        yp[j] += step
        ym[j] -= step
        fp = scaled_residual_y(yp, p)
        fm = scaled_residual_y(ym, p)
        J[:,j] = (fp - fm) / (2step)
    end
    return J
end

function refine_equilibrium(x0, p; maxiter=30, tol=1e-8)
    # protect against zeros
    x0p = max.(x0, 1e-12)
    y = log.(x0p)
    r = scaled_residual_y(y, p)
    best_y = copy(y)
    best_norm = norm(r, Inf)

    converged = false
    for iter in 1:maxiter
        r = scaled_residual_y(y, p)
        rn = norm(r, Inf)
        if rn < tol
            converged = true
            best_y = copy(y)
            best_norm = rn
            break
        end
        if rn < best_norm
            best_y = copy(y)
            best_norm = rn
        end

        Jy = jacobian_y_fd(y, p)
        # Small ridge terms handle near-singular linear systems.
        delta = try
            -(Jy + 1e-10I) \ r
        catch
            -(Jy + 1e-7I) \ r
        end

        # Limit large Newton steps.
        maxstep = maximum(abs.(delta))
        if maxstep > 2.0
            delta .*= 2.0 / maxstep
        end

        # Backtracking line search.
        accepted = false
        alpha = 1.0
        for ls in 1:12
            y_try = y + alpha * delta
            rn_try = norm(scaled_residual_y(y_try, p), Inf)
            if isfinite(rn_try) && rn_try < rn
                y = y_try
                accepted = true
                break
            end
            alpha *= 0.5
        end
        if !accepted
            # Line search stalled.
            break
        end
    end

    x_ref = exp.(best_y)
    residual_raw = norm(rhs_vec(x_ref, p), Inf)
    residual_scaled = norm(scaled_residual_y(best_y, p), Inf)
    return x_ref, converged, residual_raw, residual_scaled
end

# Up/down continuation
function run_continuation(b6_grid; direction=:up, b5=1e-4, tmax=10000.0)
    grid = direction == :up ? b6_grid : reverse(b6_grid)
    rows = DataFrame()

    # Initial condition for the first point.
    p0 = base_p(b5=b5, b6=grid[1])
    u = baseline_state(p0)

    for (k, b6) in enumerate(grid)
        println("$(direction) continuation: b6=$(b6)  ($(k)/$(length(grid)))")
        u_long, p = solve_long(b6, u; b5=b5, tmax=tmax)

        # Refine the long-time state into an equilibrium candidate.
        u_ref, conv, res_raw, res_scaled = refine_equilibrium(u_long, p)
        eigmax = max_real_eig(u_ref, p)

        push!(rows, (
            direction=String(direction),
            b5=b5,
            b6=b6,
            T=u_ref[1], M=u_ref[2], NK=u_ref[3], CTL=u_ref[4], B=u_ref[5],
            logT=log10(u_ref[1] + 1.0),
            phase=state_label(u_ref[1]),
            newton_converged=conv,
            residual_raw=res_raw,
            residual_scaled=res_scaled,
            max_real_eig=eigmax,
            stable=eigmax < 0.0
        ); cols=:union)

        # Continue from the refined state; residual columns document fit quality.
        u = u_ref
    end

    # Restore increasing b6 order for plotting and comparison.
    return sort(rows, :b6)
end

# Main run
# The range extends beyond the first hysteresis screen because the earlier
# high-tumor up-sweep branch did not collapse within the original interval.
b6_grid = logspace10(2.5e-7, 3.0e-5, 80)

println("Running equilibrium/Jacobian validation...")
println("b6 range: ", minimum(b6_grid), " to ", maximum(b6_grid))
println("grid points: ", length(b6_grid))

up_rows = run_continuation(b6_grid; direction=:up, b5=1e-4, tmax=10000.0)
down_rows = run_continuation(b6_grid; direction=:down, b5=1e-4, tmax=10000.0)
all_rows = vcat(up_rows, down_rows)
CSV.write(joinpath(outdir, "b6_equilibrium_branches.csv"), all_rows)

# Compare up/down branches at each b6.
comp = innerjoin(
    up_rows[:, [:b6, :T, :logT, :max_real_eig, :residual_scaled, :phase]],
    down_rows[:, [:b6, :T, :logT, :max_real_eig, :residual_scaled, :phase]],
    on=:b6,
    makeunique=true
)
rename!(comp, Dict(
    :T => :T_up,
    :logT => :logT_up,
    :max_real_eig => :eig_up,
    :residual_scaled => :res_up,
    :phase => :phase_up,
    :T_1 => :T_down,
    :logT_1 => :logT_down,
    :max_real_eig_1 => :eig_down,
    :residual_scaled_1 => :res_down,
    :phase_1 => :phase_down
))
comp.hysteresis_gap = abs.(comp.logT_up .- comp.logT_down)
CSV.write(joinpath(outdir, "b6_equilibrium_up_down_comparison.csv"), comp)

# Figures
# 1. Tumor branch
p1 = plot(
    up_rows.b6, safe(up_rows.T),
    xscale=:log10, yscale=:log10,
    marker=:circle, lw=2.5,
    label="up-sweep", xlabel="b6", ylabel="Tumor equilibrium candidate",
    title="Equilibrium candidates from up/down continuation",
    size=(900,600)
)
plot!(p1, down_rows.b6, safe(down_rows.T), marker=:diamond, lw=2.5, label="down-sweep")
hline!(p1, [1e3], linestyle=:dash, color=:black, label="T=1e3")
hline!(p1, [1e5], linestyle=:dot, color=:gray, label="T=1e5")
savefig(p1, joinpath(outdir, "b6_equilibrium_tumor_branches.png"))

# 2. Hysteresis gap
p2 = plot(
    comp.b6, comp.hysteresis_gap,
    xscale=:log10, marker=:circle, lw=2.5,
    xlabel="b6", ylabel="|logT_up - logT_down|",
    title="Hysteresis gap between continuation branches",
    legend=false, size=(900,600)
)
savefig(p2, joinpath(outdir, "b6_equilibrium_hysteresis_gap.png"))

# 3. Stability eigenvalue
p3 = plot(
    up_rows.b6, up_rows.max_real_eig,
    xscale=:log10, marker=:circle, lw=2.5,
    xlabel="b6", ylabel="max real eigenvalue",
    title="Local stability screen at equilibrium candidates",
    label="up-sweep", size=(900,600)
)
plot!(p3, down_rows.b6, down_rows.max_real_eig, marker=:diamond, lw=2.5, label="down-sweep")
hline!(p3, [0.0], linestyle=:dash, color=:black, label="0")
savefig(p3, joinpath(outdir, "b6_equilibrium_stability.png"))

# 4. Residuals
p4 = plot(
    up_rows.b6, safe(up_rows.residual_scaled),
    xscale=:log10, yscale=:log10,
    marker=:circle, lw=2.5,
    xlabel="b6", ylabel="scaled residual ||f(x*)||",
    title="Steady-state residual check",
    label="up-sweep", size=(900,600)
)
plot!(p4, down_rows.b6, safe(down_rows.residual_scaled), marker=:diamond, lw=2.5, label="down-sweep")
hline!(p4, [1e-6, 1e-8], linestyle=:dash, color=[:gray :black], label=["1e-6" "1e-8"])
savefig(p4, joinpath(outdir, "b6_equilibrium_residuals.png"))

println("Done.")
println("Outputs saved to: ", outdir)
println("Key files:")
println("  b6_equilibrium_branches.csv")
println("  b6_equilibrium_up_down_comparison.csv")
println("  b6_equilibrium_tumor_branches.png")
println("  b6_equilibrium_hysteresis_gap.png")
println("  b6_equilibrium_stability.png")
println("  b6_equilibrium_residuals.png")
