using DifferentialEquations
using Plots
using Colors
using Statistics
using DataFrames
using CSV

logo = Colors.JULIA_LOGO_COLORS
cell_colors = [logo.red, colorant"goldenrod2", logo.green, logo.blue, colorant"purple"]

# Compare the two CTL formulations near the b6 threshold

# direct modeL: B_to_CTL = a11 * xB / (g5 + xB)
# coupled mode: B_to_CTL = a11 * (xT^2 / (g3 + xT^2)) * (xB / (g5 + xB))

function base_p(; a10 = 2.0, b5 = 1.0e-4, a11 = 1.0e-1, b6 = 1.0e-7,
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

# start near the tumor-free equilibrium
function tumor_free_start(p; tumor_cells = 2.0)
    xM0 = p.a2 / p.z2
    xN0 = p.z2 * p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC0 = 0.0
    xB0 = p.a8 / (p.z5 + p.b5 * xM0)
    return [tumor_cells, xM0, xN0, xC0, xB0]
end

function model_B!(du, u, p, t)
    xT, xM, xN, xC, xB = u
    # switch between the two CTL help formulations
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

function solve_ode(; a10 = 2.0, b5 = 1e-4, a11 = 1e-1, b6 = 1e-7,
                   b1 = 3.5e-6, b2 = 1.1e-7, ctl_mode = :direct,
                   tmax = 365.0, saveat = 1.0)
    p = base_p(a10 = a10, b5 = b5, a11 = a11, b6 = b6, b1 = b1, b2 = b2, ctl_mode = ctl_mode)
    u0 = tumor_free_start(p; tumor_cells = 2.0)
    prob = ODEProblem(model_B!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5(); saveat = saveat, reltol = 1e-7, abstol = 1e-9,
                isoutofdomain = (u,p,t)->any(x->x<0,u))
    return sol
end

function trapz(x, y)
    # trapezoid rule for AUC
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
    #Summary numbers from one trajectory
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

function safe_log_vec(v; eps=1e-6)
    max.(collect(v), eps)
end

function compare_ctl_forms(; a10=2.0, b5=1e-4, a11=1e-1)
    b6_vals = [2e-7, 3e-7, 4e-7, 5e-7, 7e-7, 1e-6, 1.5e-6, 2e-6]
    modes = [:direct, :coupled]

    df = DataFrame(mode=String[], b6=Float64[], T60=Float64[], T100=Float64[], T365=Float64[],
                   M365=Float64[], NK365=Float64[], CTL365=Float64[], B365=Float64[],
                   t_to_1e5=Float64[], t_to_1e6=Float64[], auc_0_120=Float64[])

    for mode in modes
        for b6 in b6_vals
            sol = solve_ode(a10=a10, b5=b5, a11=a11, b6=b6, ctl_mode=mode)
            m = metrics_from_sol(sol)
            push!(df, (String(mode), b6, m.T60, m.T100, m.T365, m.M365, m.NK365, m.CTL365, m.B365, m.t_to_1e5, m.t_to_1e6, m.auc_0_120))
        end
    end
    CSV.write("ctl_form_compare_near_b6.csv", df)

    # summary plots
    function mode_df(mode)
        filter(row -> row.mode == String(mode), df)
    end
    ddir = mode_df(:direct)
    dcpl = mode_df(:coupled)

    p1 = plot(ddir.b6, ddir.T365, xscale=:log10, lw=2.5, marker=:circle, label="direct", xlabel="b6", ylabel="Tumor day 365", title="Tumor day 365")
    plot!(p1, dcpl.b6, dcpl.T365, xscale=:log10, lw=2.5, marker=:diamond, label="coupled")

    p2 = plot(ddir.b6, ddir.T100, xscale=:log10, lw=2.5, marker=:circle, label="direct", xlabel="b6", ylabel="Tumor day 100", title="Tumor day 100")
    plot!(p2, dcpl.b6, dcpl.T100, xscale=:log10, lw=2.5, marker=:diamond, label="coupled")

    p3 = plot(ddir.b6, ddir.CTL365, xscale=:log10, lw=2.5, marker=:circle, label="direct", xlabel="b6", ylabel="CTL day 365", title="CTL day 365")
    plot!(p3, dcpl.b6, dcpl.CTL365, xscale=:log10, lw=2.5, marker=:diamond, label="coupled")

    p4 = plot(ddir.b6, ddir.t_to_1e5, xscale=:log10, lw=2.5, marker=:circle, label="direct", xlabel="b6", ylabel="Time to 1e5", title="Time to 1e5")
    plot!(p4, dcpl.b6, dcpl.t_to_1e5, xscale=:log10, lw=2.5, marker=:diamond, label="coupled")

    summary = plot(p1,p2,p3,p4, layout=(2,2), size=(1100,850))
    savefig(summary, "ctl_form_compare_near_b6_summary.png")

    # trajectory panels
    for mode in modes
        subplots = Any[]
        for b6 in [4e-7, 5e-7, 7e-7, 1e-6]
            sol = solve_ode(a10=a10, b5=b5, a11=a11, b6=b6, ctl_mode=mode)
            p = plot(sol.t, safe_log_vec(sol[1,:]), yscale=:log10, color=cell_colors[1], lw=2.2,
                     xlabel="Time (days)", ylabel="Cells", label="Tumor",
                     title="mode=$(mode), b6=$(b6)")
            plot!(p, sol.t, safe_log_vec(sol[2,:]), yscale=:log10, color=cell_colors[2], lw=2.0, label="MDSC")
            plot!(p, sol.t, safe_log_vec(sol[3,:]), yscale=:log10, color=cell_colors[3], lw=2.0, label="NK")
            plot!(p, sol.t, safe_log_vec(sol[4,:]), yscale=:log10, color=cell_colors[4], lw=2.0, label="CTL")
            plot!(p, sol.t, safe_log_vec(sol[5,:]), yscale=:log10, color=cell_colors[5], lw=2.0, label="B")
            push!(subplots, p)
        end
        fig = plot(subplots..., layout=(length(subplots),1), size=(1000, 250*length(subplots)))
        savefig(fig, "ctl_form_compare_$(mode)_trajectories.png")
    end

    println("Saved ctl_form_compare_near_b6.csv")
    println("Saved ctl_form_compare_near_b6_summary.png")
    println("Saved ctl_form_compare_direct_trajectories.png")
    println("Saved ctl_form_compare_coupled_trajectories.png")

    return df
end

# run
compare_ctl_forms(a10=2.0, b5=1e-4, a11=1e-1)
println("SUCESS!")
