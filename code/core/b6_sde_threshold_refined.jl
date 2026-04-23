using DifferentialEquations
using Random
using Statistics
using DataFrames
using CSV
using Plots

# ------------------------------------------------------------
# Focused SDE refinement for b6 near the transition region.
# Goal: get stochastic results that are easy to explain.
#
# Main changes vs the earlier SDE script:
# 1) keep the simpler DIRECT CTL formulation
# 2) remove the terminate-on-negative callback
# 3) instead, clamp states to zero after each step
# 4) summarize threshold-crossing probabilities instead of the old
#    "survives above 1 cell" rule
# ------------------------------------------------------------

function base_p(; a10 = 2.0, b5 = 1.0e-4, a11 = 1.0e-1, b6 = 5.0e-7,
                 b1 = 3.5e-6, b2 = 1.1e-7, ctl_mode = :direct)
    return (
        a1 = 1.0e-1,
        h  = 1.0e7,
        b1 = b1,
        b2 = b2,
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
        b5 = b5,
        z5 = 2.0e-2,
        g4 = 2.02e7,
        g5 = 1.0e3,
        a10 = a10,
        a11 = a11,
        b6 = b6,
        ctl_mode = ctl_mode,
    )
end

function tumor_free_start(p; tumor_cells = 2.0)
    xM0 = p.a2 / p.z2
    xN0 = p.z2 * p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC0 = 0.0
    xB0 = p.a8 / (p.z5 + p.b5 * xM0)
    return [tumor_cells, xM0, xN0, xC0, xB0]
end

function model_B!(du, u, p, t)
    xT, xM, xN, xC, xB = max.(u, 0.0)

    B_to_CTL = if p.ctl_mode == :direct
        p.a11 * (xB / (p.g5 + xB))
    else
        p.a11 * (xT^2 / (p.g3 + xT^2)) * (xB / (p.g5 + xB))
    end

    du[1] = p.a1 * xT * log(max(p.h / max(xT, 1e-12), 1.0)) - p.b1 * xT * xN - p.b2 * xT * xC - p.b6 * xT * xB - p.z1 * xT
    du[2] = p.a2 + p.a3 * xT / (p.g1 + xT) - p.z2 * xM
    du[3] = p.a4 + p.a5 * (xT^2) / (p.g2 + xT^2) * (1.0 + p.a10 * xB / (p.g5 + xB)) - p.b3 * xM * xN - p.z3 * xN
    du[4] = p.a6 * xT * xN + p.a7 * (xT^2) / (p.g3 + xT^2) + B_to_CTL - p.b4 * xM * xC - p.z4 * xC
    du[5] = p.a8 + p.a9 * (xT^2) / (p.g4 + xT^2) - p.b5 * xM * xB - p.z5 * xB
end

function sigma_B!(du, u, p, t)
    for i in eachindex(u)
        du[i] = max(u[i], 0.0)
    end
end

function clamp_negative!(integrator)
    @inbounds for i in eachindex(integrator.u)
        if integrator.u[i] < 0.0
            integrator.u[i] = 0.0
        end
    end
end
cb_clamp = FunctionCallingCallback(clamp_negative!; func_everystep = true)

function solve_sde_once(; a10 = 2.0, b5 = 1e-4, a11 = 1e-1, b6 = 5e-7,
                        b1 = 3.5e-6, b2 = 1.1e-7, ctl_mode = :direct,
                        tmax = 365.0, dt = 0.02, saveat = 1.0, rngseed = 1)
    Random.seed!(rngseed)
    p = base_p(a10 = a10, b5 = b5, a11 = a11, b6 = b6, b1 = b1, b2 = b2, ctl_mode = ctl_mode)
    u0 = tumor_free_start(p; tumor_cells = 2.0)
    prob = SDEProblem(model_B!, sigma_B!, u0, (0.0, tmax), p)
    sol = solve(prob, SOSRI(); dt = dt, saveat = saveat,
                reltol = 1e-5, abstol = 1e-7, callback = cb_clamp)
    return sol
end

function first_passage_time(v, t, thresh)
    idx = findfirst(x -> x >= thresh, v)
    return isnothing(idx) ? NaN : t[idx]
end

function summarize_run(sol)
    T = sol[1,:]
    t = sol.t
    return (
        Tfinal = T[end],
        Tmax = maximum(T),
        ptime_1e2 = first_passage_time(T, t, 1e2),
        ptime_1e3 = first_passage_time(T, t, 1e3),
        extinct = T[end] < 1.0,
        reach_1e2 = maximum(T) >= 1e2,
        reach_1e3 = maximum(T) >= 1e3,
        final_ge_1e2 = T[end] >= 1e2,
        final_ge_1e3 = T[end] >= 1e3,
    )
end

function wilson_interval(k, n; z=1.96)
    p = k / n
    denom = 1 + z^2/n
    center = (p + z^2/(2n)) / denom
    half = z * sqrt((p*(1-p) + z^2/(4n))/n) / denom
    return (center-half, center+half)
end

function run_b6_sde_refined(; a10 = 2.0, a11 = 1e-1, ctl_mode = :direct,
                             b5_vals = [1e-5, 1e-4, 1e-3],
                             b6_vals = [4e-7, 5e-7, 6e-7, 7e-7, 8e-7, 1e-6],
                             nsims = 600,
                             dt = 0.02)
    rows = DataFrame(
        b5 = Float64[], b6 = Float64[], nsims = Int[],
        p_reach_1e2 = Float64[], p_reach_1e3 = Float64[],
        p_final_ge_1e2 = Float64[], p_final_ge_1e3 = Float64[],
        p_extinct = Float64[],
        p_reach_1e2_lo = Float64[], p_reach_1e2_hi = Float64[],
        p_final_ge_1e2_lo = Float64[], p_final_ge_1e2_hi = Float64[],
        mean_Tfinal = Float64[], median_Tfinal = Float64[],
        mean_Tmax = Float64[], median_Tmax = Float64[],
        mean_t_1e2 = Float64[], mean_t_1e3 = Float64[]
    )

    for b5 in b5_vals
        for b6 in b6_vals
            finals = Float64[]
            tmaxes = Float64[]
            t1e2 = Float64[]
            t1e3 = Float64[]
            n_reach_1e2 = 0
            n_reach_1e3 = 0
            n_final_ge_1e2 = 0
            n_final_ge_1e3 = 0
            n_extinct = 0

            for s in 1:nsims
                sol = solve_sde_once(a10=a10, b5=b5, a11=a11, b6=b6, ctl_mode=ctl_mode, rngseed=s, dt=dt)
                m = summarize_run(sol)
                push!(finals, m.Tfinal)
                push!(tmaxes, m.Tmax)
                if !isnan(m.ptime_1e2)
                    push!(t1e2, m.ptime_1e2)
                    n_reach_1e2 += 1
                end
                if !isnan(m.ptime_1e3)
                    push!(t1e3, m.ptime_1e3)
                    n_reach_1e3 += 1
                end
                n_final_ge_1e2 += m.final_ge_1e2 ? 1 : 0
                n_final_ge_1e3 += m.final_ge_1e3 ? 1 : 0
                n_extinct += m.extinct ? 1 : 0
            end

            lo1, hi1 = wilson_interval(n_reach_1e2, nsims)
            lo2, hi2 = wilson_interval(n_final_ge_1e2, nsims)

            push!(rows, (
                b5, b6, nsims,
                n_reach_1e2/nsims, n_reach_1e3/nsims,
                n_final_ge_1e2/nsims, n_final_ge_1e3/nsims,
                n_extinct/nsims,
                lo1, hi1, lo2, hi2,
                mean(finals), median(finals), mean(tmaxes), median(tmaxes),
                isempty(t1e2) ? NaN : mean(t1e2),
                isempty(t1e3) ? NaN : mean(t1e3)
            ))
        end
    end

    CSV.write("b6_sde_threshold_refined.csv", rows)

    p1 = plot(xscale=:log10, xlabel="b6", ylabel="P(max T ≥ 1e2)", title="Chance of meaningful outgrowth")
    p2 = plot(xscale=:log10, xlabel="b6", ylabel="P(final T ≥ 1e2)", title="Chance tumor still established at day 365")
    p3 = plot(xscale=:log10, xlabel="b6", ylabel="Median final tumor", title="Median final tumor size")
    p4 = plot(xscale=:log10, xlabel="b6", ylabel="Median max tumor", title="Median max tumor size")
    p5 = plot(xscale=:log10, xlabel="b6", ylabel="Mean time to 1e2", title="Mean time to reach 1e2")
    p6 = plot(xscale=:log10, xlabel="b6", ylabel="P(extinct by day 365)", title="Extinction probability")

    for b5 in b5_vals
        sub = rows[rows.b5 .== b5, :]
        lab = "b5=$(b5)"
        plot!(p1, sub.b6, sub.p_reach_1e2, marker=:circle, lw=2, label=lab)
        plot!(p2, sub.b6, sub.p_final_ge_1e2, marker=:circle, lw=2, label=lab)
        plot!(p3, sub.b6, sub.median_Tfinal, marker=:circle, lw=2, label=lab)
        plot!(p4, sub.b6, sub.median_Tmax, marker=:circle, lw=2, label=lab)
        plot!(p5, sub.b6, sub.mean_t_1e2, marker=:circle, lw=2, label=lab)
        plot!(p6, sub.b6, sub.p_extinct, marker=:circle, lw=2, label=lab)
    end

    fig = plot(p1,p2,p3,p4,p5,p6, layout=(3,2), size=(1100,1200))
    savefig(fig, "b6_sde_threshold_refined.png")

    return rows
end

println("Running refined SDE around the b6 threshold...")
rows = run_b6_sde_refined()
println("Saved b6_sde_threshold_refined.csv and b6_sde_threshold_refined.png")
