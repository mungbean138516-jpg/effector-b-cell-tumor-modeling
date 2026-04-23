using DifferentialEquations
using Plots
using Colors
using Statistics
using DataFrames
using CSV
using Random

# beta 6 analysis for the effector-B model

# baseline parameters

function base_p(; a10 = 2.0, b5 = 1.0e-4, a11 = 1.0e-1, b6 = 1.0e-7,
                 b1 = 3.5e-6, b2 = 1.1e-7, ctl_mode = :direct)
    return (
        # original CIR parameters
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

        # B-cell extension
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

# tumor-free immune baseline
function tumor_free_start(p; tumor_cells = 2.0)
    xM0 = p.a2 / p.z2
    xN0 = p.z2 * p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC0 = 0.0
    xB0 = p.a8 / (p.z5 + p.b5 * xM0)
    return [tumor_cells, xM0, xN0, xC0, xB0]
end

# ODE+SDE system
# 
function model_B!(du, u, p, t)
    xT, xM, xN, xC, xB = u

    B_to_CTL = if p.ctl_mode == :direct
        p.a11 * (xB / (p.g5 + xB))
    else
        p.a11 * (xT^2 / (p.g3 + xT^2)) * (xB / (p.g5 + xB))
    end

    # tumor
    du[1] = p.a1 * xT * log(max(p.h / max(xT, 1e-12), 1.0)) - p.b1 * xT * xN - p.b2 * xT * xC - p.b6 * xT * xB - p.z1 * xT

    # MDSC
    du[2] = p.a2 + p.a3 * xT / (p.g1 + xT) - p.z2 * xM

    # NK
    du[3] = p.a4 + p.a5 * (xT^2) / (p.g2 + xT^2) * (1.0 + p.a10 * xB / (p.g5 + xB)) - p.b3 * xM * xN - p.z3 * xN

    # CTL
    du[4] = p.a6 * xT * xN + p.a7 * (xT^2) / (p.g3 + xT^2) + B_to_CTL - p.b4 * xM * xC - p.z4 * xC

    # B
    du[5] = p.a8 + p.a9 * (xT^2) / (p.g4 + xT^2) - p.b5 * xM * xB - p.z5 * xB
end

function sigma_B!(du, u, p, t)
    for i in eachindex(u)
        # multiplicative noise
        du[i] = max(u[i], 0.0)  
    end
end

function stop_negative(u, t, integrator)
    any(x -> x < 0.0, u)
end
cb_stop_negative = DiscreteCallback(stop_negative, terminate!)

function solve_ode(; a10 = 2.0, b5 = 1e-4, a11 = 1e-1, b6 = 1e-7,
                   b1 = 3.5e-6, b2 = 1.1e-7, ctl_mode = :direct,
                   tmax = 365.0, saveat = 1.0)
    p = base_p(a10 = a10, b5 = b5, a11 = a11, b6 = b6, b1 = b1, b2 = b2, ctl_mode = ctl_mode)
    u0 = tumor_free_start(p; tumor_cells = 2.0)
    prob = ODEProblem(model_B!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5(); saveat = saveat, reltol = 1e-7, abstol = 1e-9,
                callback = cb_stop_negative,
                isoutofdomain = (u,p,t)->any(x->x<0,u))
    return sol, p
end

function solve_sde_once(; a10 = 2.0, b5 = 1e-4, a11 = 1e-1, b6 = 1e-7,
                        b1 = 3.5e-6, b2 = 1.1e-7, ctl_mode = :direct,
                        tmax = 365.0, dt = 0.1, rngseed = 1)
    Random.seed!(rngseed)
    p = base_p(a10 = a10, b5 = b5, a11 = a11, b6 = b6, b1 = b1, b2 = b2, ctl_mode = ctl_mode)
    u0 = tumor_free_start(p; tumor_cells = 2.0)
    prob = SDEProblem(model_B!, sigma_B!, u0, (0.0, tmax), p)
    sol = solve(prob, SOSRI(); dt = dt, saveat = 1.0,
                reltol = 1e-5, abstol = 1e-7,
                callback = cb_stop_negative,
                isoutofdomain = (u,p,t)->any(x->x<0,u))
    return sol, p
end

# helpers
function trapz(x, y)
    s = 0.0
    for i in 1:length(x)-1
        s += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
    end
    return s
end

function time_to_threshold(sol, thresh)
    idx = findfirst(x -> x >= thresh, sol[1,:])
    return isnothing(idx) ? NaN : sol.t[idx]
end

function metrics_from_sol(sol)
    t = sol.t
    xT = sol[1,:]
    xM = sol[2,:]
    xN = sol[3,:]
    xC = sol[4,:]
    xB = sol[5,:]

    idx60 = findfirst(==(60.0), t)
    idx100 = findfirst(==(100.0), t)

    return (
        T60 = isnothing(idx60) ? NaN : xT[idx60],
        T100 = isnothing(idx100) ? NaN : xT[idx100],
        T365 = xT[end],
        M365 = xM[end],
        NK365 = xN[end],
        CTL365 = xC[end],
        B365 = xB[end],
        t_to_1e5 = time_to_threshold(sol, 1e5),
        t_to_1e6 = time_to_threshold(sol, 1e6),
        auc_0_120 = trapz(t[1:121], xT[1:121])
    )
end

function safe_log_vec(v; eps = 1e-6)
    return max.(collect(v), eps)
end

# ODE b6 sweep
function run_refined_b6_sweep(; a10 = 2.0, b5 = 1e-4, a11 = 1e-1, ctl_mode=:direct)
    b6_vals = [1e-7, 2e-7, 3e-7, 4e-7, 5e-7, 7e-7, 1e-6, 1.5e-6, 2e-6, 3e-6, 5e-6]

    df = DataFrame(b6 = Float64[], T60=Float64[], T100=Float64[], T365=Float64[],
                   M365=Float64[], NK365=Float64[], CTL365=Float64[], B365=Float64[],
                   t_to_1e5=Float64[], t_to_1e6=Float64[], auc_0_120=Float64[])

    plt = plot(layout=(length(b6_vals),1), size=(1000, 250*length(b6_vals)))

    for (k, b6) in enumerate(b6_vals)
        sol, _ = solve_ode(a10=a10, b5=b5, a11=a11, b6=b6, ctl_mode=ctl_mode)
        m = metrics_from_sol(sol)
        push!(df, (b6, m.T60, m.T100, m.T365, m.M365, m.NK365, m.CTL365, m.B365, m.t_to_1e5, m.t_to_1e6, m.auc_0_120))

        psub = plot(sol.t, safe_log_vec(sol[1,:]), yscale=:log10, color=cell_colors[1], lw=2.2,
                    xlabel="Time (days)", ylabel="Cells", label="Tumor",
                    title="ODE trajectories, b6 = $(b6), b5 = $(b5)")
        plot!(psub, sol.t, safe_log_vec(sol[2,:]), yscale=:log10, color=cell_colors[2], lw=2.0, label="MDSC")
        plot!(psub, sol.t, safe_log_vec(sol[3,:]), yscale=:log10, color=cell_colors[3], lw=2.0, label="NK")
        plot!(psub, sol.t, safe_log_vec(sol[4,:]), yscale=:log10, color=cell_colors[4], lw=2.0, label="CTL")
        plot!(psub, sol.t, safe_log_vec(sol[5,:]), yscale=:log10, color=cell_colors[5], lw=2.0, label="B")
        plot!(plt[k], psub)
    end

    savefig(plt, "refined_b6_ode_trajectories.png")
    CSV.write("refined_b6_ode_sweep.csv", df)

    p1 = plot(df.b6, df.T100, xscale=:log10, marker=:circle, lw=2.5, xlabel="b6", ylabel="Tumor at day 100", title="Tumor day 100 vs b6")
    p2 = plot(df.b6, df.T365, xscale=:log10, marker=:circle, lw=2.5, xlabel="b6", ylabel="Tumor at day 365", title="Tumor day 365 vs b6")
    p3 = plot(df.b6, df.t_to_1e5, xscale=:log10, marker=:circle, lw=2.5, xlabel="b6", ylabel="Time to 1e5", title="Time to 1e5 vs b6")
    p4 = plot(df.b6, df.auc_0_120, xscale=:log10, marker=:circle, lw=2.5, xlabel="b6", ylabel="Tumor AUC 0-120 d", title="Early tumor burden vs b6")
    p5 = plot(df.b6, df.B365, xscale=:log10, marker=:circle, lw=2.5, xlabel="b6", ylabel="B at day 365", title="B abundance vs b6")
    p6 = plot(df.b6, df.NK365, xscale=:log10, marker=:circle, lw=2.5, xlabel="b6", ylabel="NK at day 365", title="NK abundance vs b6")
    summ = plot(p1,p2,p3,p4,p5,p6, layout=(3,2), size=(1100,1200))
    savefig(summ, "refined_b6_ode_summary.png")

    return df
end

# b5*b6 heatmaps
function make_b5_b6_heatmaps(; a10 = 2.0, a11 = 1e-1, ctl_mode=:direct)
    b5_vals = [1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
    b6_vals = [1e-7, 2e-7, 3e-7, 4e-7, 5e-7, 7e-7, 1e-6, 1.5e-6, 2e-6, 3e-6, 5e-6]

    Z100 = zeros(length(b5_vals), length(b6_vals))
    Z365 = zeros(length(b5_vals), length(b6_vals))
    Zt = fill(NaN, length(b5_vals), length(b6_vals))

    for (i, b5) in enumerate(b5_vals)
        for (j, b6) in enumerate(b6_vals)
            sol, _ = solve_ode(a10=a10, b5=b5, a11=a11, b6=b6, ctl_mode=ctl_mode)
            m = metrics_from_sol(sol)
            Z100[i,j] = m.T100
            Z365[i,j] = m.T365
            Zt[i,j] = isnan(m.t_to_1e5) ? 400.0 : m.t_to_1e5
        end
    end

    p100 = heatmap(b6_vals, b5_vals, Z100, xscale=:log10, yscale=:log10,
                   xlabel="b6", ylabel="b5", title="Tumor at day 100 (b5 x b6)", colorbar_title="T100")
    savefig(p100, "heatmap_T100_b5_b6.png")

    p365 = heatmap(b6_vals, b5_vals, Z365, xscale=:log10, yscale=:log10,
                   xlabel="b6", ylabel="b5", title="Tumor at day 365 (b5 x b6)", colorbar_title="T365")
    savefig(p365, "heatmap_T365_b5_b6.png")

    pt = heatmap(b6_vals, b5_vals, Zt, xscale=:log10, yscale=:log10,
                 xlabel="b6", ylabel="b5", title="Time to 1e5 tumor cells (b5 x b6)", colorbar_title="days")
    savefig(pt, "heatmap_t1e5_b5_b6.png")

    return (b5_vals=b5_vals, b6_vals=b6_vals, Z100=Z100, Z365=Z365, Zt=Zt)
end

# SDE
function run_refined_b6_sde(; a10 = 2.0, a11 = 1e-1, ctl_mode=:direct,
                            b5_vals = [1e-5, 1e-4, 1e-3],
                            b6_vals = [2e-7, 5e-7, 1e-6, 2e-6],
                            nsims = 300)
    rows = DataFrame(b5=Float64[], b6=Float64[], success_prob=Float64[], mean_success_tumor=Float64[], mean_time_to_extinction=Float64[])

    for b5 in b5_vals
        for b6 in b6_vals
            successes = Bool[]
            success_sizes = Float64[]
            ext_times = Float64[]

            for s in 1:nsims
                sol, _ = solve_sde_once(a10=a10, b5=b5, a11=a11, b6=b6, ctl_mode=ctl_mode, rngseed=s)
                tumor = sol[1,:]
                ok = minimum(floor.(tumor)) > 0.0
                push!(successes, ok)
                if ok
                    push!(success_sizes, mean(tumor))
                else
                    idx = findfirst(x -> floor(x) < 1.0, tumor)
                    if isnothing(idx)
                        push!(ext_times, sol.t[end])
                    else
                        push!(ext_times, sol.t[idx])
                    end
                end
            end

            push!(rows, (b5, b6,
                         mean(successes),
                         isempty(success_sizes) ? NaN : mean(success_sizes),
                         isempty(ext_times) ? NaN : mean(ext_times)))
        end
    end

    CSV.write("refined_b6_sde_summary.csv", rows)

    #quick plots by b6 for each b5
    plt1 = plot(xscale=:log10, xlabel="b6", ylabel="Success probability", title="SDE metastasis establishment vs b6")
    plt2 = plot(xscale=:log10, xlabel="b6", ylabel="Mean tumor size (successful runs)", title="Successful-run tumor burden vs b6")
    plt3 = plot(xscale=:log10, xlabel="b6", ylabel="Mean time to extinction", title="Extinction time vs b6")

    for b5 in b5_vals
        sub = rows[rows.b5 .== b5, :]
        lab = "b5=$(b5)"
        plot!(plt1, sub.b6, sub.success_prob, marker=:circle, lw=2.2, label=lab)
        plot!(plt2, sub.b6, sub.mean_success_tumor, marker=:circle, lw=2.2, label=lab)
        plot!(plt3, sub.b6, sub.mean_time_to_extinction, marker=:circle, lw=2.2, label=lab)
    end

    savefig(plot(plt1, plt2, plt3, layout=(3,1), size=(900,1200)), "refined_b6_sde_summary.png")
    return rows
end

# main
println("Running refined b6 ODE sweep...")
odf = run_refined_b6_sweep()
println("Saved refined_b6_ode_sweep.csv and ODE figures.")

println("Making b5 x b6 heatmaps...")
heat = make_b5_b6_heatmaps()
println("Saved b5 x b6 heatmaps.")

println("Running focused SDE checks around the b6 threshold...")
sdf = run_refined_b6_sde()
println("Saved refined_b6_sde_summary.csv and SDE figure.")

println("Done.!")
