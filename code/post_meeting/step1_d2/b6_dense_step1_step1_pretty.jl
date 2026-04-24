using DifferentialEquations
using Plots
using DataFrames
using CSV
using Statistics

# finer grade, 10000 days

default(
    legendfontsize = 9,
    guidefontsize = 11,
    tickfontsize = 9,
    titlefontsize = 14,
    linewidth = 2.4,
    markersize = 4.5,
    framestyle = :box,
    grid = true,
    gridalpha = 0.25,
    minorgrid = false,
)

outdir = "b6_dense_step1_outputs_pretty"
mkpath(outdir)

# Small helpers

function logspace10(lo, hi, n)
    return 10 .^ range(log10(lo), log10(hi), length=n)
end

safe(v; eps=1e-8) = max.(collect(v), eps)

function sci_label(x)
    return string(round(x, sigdigits=3))
end

function grid_tag(x, coarse_vals, fine_vals)
    in_coarse = any(isapprox(x, y; rtol=1e-12, atol=0.0) for y in coarse_vals)
    in_fine   = any(isapprox(x, y; rtol=1e-12, atol=0.0) for y in fine_vals)
    if in_coarse && in_fine
        return "both"
    elseif in_fine
        return "fine"
    else
        return "coarse"
    end
end

function nearest_value(vals, target)
    return vals[argmin(abs.(vals .- target))]
end

# Fixed baseline model

function base_p(; b6=1e-7)
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

        # effector-B extension
        a8  = 1.0e4,
        a9  = 5.0e-2,
        b5  = 1.0e-4, # fixed in Step 1
        z5  = 2.0e-2,
        g4  = 2.02e7,
        g5  = 1.0e3,
        a10 = 2.0,
        a11 = 1.0e-1,
        b6  = b6,
    )
end

function start_state(p; tumor_cells=2.0)
    xM0 = p.a2 / p.z2
    xN0 = p.z2 * p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC0 = 0.0
    xB0 = p.a8 / (p.z5 + p.b5 * xM0)
    return [tumor_cells, xM0, xN0, xC0, xB0]
end

# direct CTL
function model!(du, u, p, t)
    xT, xM, xN, xC, xB = max.(u, 0.0)

    B_to_CTL = p.a11 * xB / (p.g5 + xB)

    growth = p.a1 * xT * log(max(p.h / max(xT, 1e-12), 1.0))
    du[1] = growth - p.b1*xT*xN - p.b2*xT*xC - p.b6*xT*xB - p.z1*xT
    du[2] = p.a2 + p.a3*xT/(p.g1 + xT) - p.z2*xM

    nk_recruit = p.a5 * xT^2/(p.g2 + xT^2) * (1.0 + p.a10*xB/(p.g5 + xB))
    du[3] = p.a4 + nk_recruit - p.b3*xM*xN - p.z3*xN

    du[4] = p.a6*xT*xN + p.a7*xT^2/(p.g3 + xT^2) + B_to_CTL - p.b4*xM*xC - p.z4*xC
    du[5] = p.a8 + p.a9*xT^2/(p.g4 + xT^2) - p.b5*xM*xB - p.z5*xB
end

function bad_state(u, t, integrator)
    any(x -> !isfinite(x) || x < -1e-6, u)
end
cb_bad = DiscreteCallback(bad_state, terminate!)

function solve_one_b6(b6; tmax=10000.0, saveat=1.0)
    p = base_p(b6=b6)
    u0 = start_state(p)
    prob = ODEProblem(model!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5();
                saveat=saveat,
                reltol=1e-7,
                abstol=1e-9,
                callback=cb_bad,
                isoutofdomain=(u,p,t)->any(x->x < -1e-6, u))
    return sol
end

# Metrics

function nearest_index(times, day)
    return argmin(abs.(times .- day))
end

function value_at(sol, state_index, day)
    idx = nearest_index(sol.t, day)
    return sol[state_index, idx]
end

function trapz(times, vals)
    s = 0.0
    for i in 1:(length(times)-1)
        s += 0.5 * (vals[i] + vals[i+1]) * (times[i+1] - times[i])
    end
    return s
end

function auc_until(sol, day)
    idxs = findall(t -> t <= day, sol.t)
    return trapz(sol.t[idxs], sol[1,idxs])
end

function time_to_threshold(sol, thresh)
    idx = findfirst(x -> x >= thresh, sol[1,:])
    return isnothing(idx) ? NaN : sol.t[idx]
end

function state_label(Tfinal)
    if Tfinal < 1e3
        return "tumor_controlled"
    elseif Tfinal > 1e5
        return "tumor_dominant"
    else
        return "intermediate"
    end
end

function metrics(sol, b6, coarse_vals, fine_vals)
    T100    = value_at(sol, 1, 100.0)
    T365    = value_at(sol, 1, 365.0)
    T10000  = value_at(sol, 1, 10000.0)
    M10000  = value_at(sol, 2, 10000.0)
    NK10000 = value_at(sol, 3, 10000.0)
    CTL10000= value_at(sol, 4, 10000.0)
    B10000  = value_at(sol, 5, 10000.0)

    t1e5 = time_to_threshold(sol, 1e5)
    t1e6 = time_to_threshold(sol, 1e6)

    return (
        b6=b6,
        grid_type=grid_tag(b6, coarse_vals, fine_vals),
        T100=T100,
        T365=T365,
        T10000=T10000,
        M10000=M10000,
        NK10000=NK10000,
        CTL10000=CTL10000,
        B10000=B10000,
        reached_1e5=!isnan(t1e5),
        reached_1e6=!isnan(t1e6),
        t_to_1e5=t1e5,
        t_to_1e6=t1e6,
        t_to_1e5_plot=isnan(t1e5) ? 10000.0 : t1e5,
        t_to_1e6_plot=isnan(t1e6) ? 10000.0 : t1e6,
        auc_0_365=auc_until(sol, 365.0),
        auc_0_1000=auc_until(sol, 1000.0),
        final_state=state_label(T10000)
    )
end


# Dense b6 grids

coarse_b6 = logspace10(1e-7, 3e-6, 40)
fine_lo = 4e-7
fine_hi = 9e-7
fine_b6   = collect(range(fine_lo, fine_hi, length=40))
all_b6    = sort(unique(vcat(coarse_b6, fine_b6)))

# Representative values for the main trajectories figure
rep_global_targets = [1e-7, 2e-7, 4e-7, 5e-7, 6e-7, 7e-7, 8e-7, 1e-6, 2e-6, 3e-6]
rep_global = unique([nearest_value(all_b6, x) for x in rep_global_targets])


rep_thresh_targets = [4.0e-7, 4.5e-7, 5.0e-7, 5.5e-7, 6.0e-7, 6.5e-7, 7.0e-7, 7.5e-7, 8.0e-7, 8.5e-7, 9.0e-7]
rep_thresh = unique([nearest_value(all_b6, x) for x in rep_thresh_targets])

println("Step 0 fixed:")
println("  direct CTL formulation")
println("  b5 fixed at 1e-4")
println("  baseline initial condition")
println("  varying b6 only")
println("Step 1 dense sweep:")
println("  coarse points = ", length(coarse_b6))
println("  fine points   = ", length(fine_b6))
println("  total unique  = ", length(all_b6))

# Run the sweep

rows = DataFrame()
solutions_for_reps = Dict{Float64, Any}()

for (k, b6) in enumerate(all_b6)
    println("Running b6 = ", b6, "  (", k, "/", length(all_b6), ")")
    sol = solve_one_b6(b6)
    push!(rows, metrics(sol, b6, coarse_b6, fine_b6); cols=:union)
end
sort!(rows, :b6)
CSV.write(joinpath(outdir, "b6_dense_step1_sweep_pretty.csv"), rows)
println("Saved CSV.")

# only the representative runs for plotting
for b6 in unique(vcat(rep_global, rep_thresh))
    solutions_for_reps[b6] = solve_one_b6(b6)
end

# Plot settings
cell_names = ["Tumor", "MDSC", "NK", "CTL", "B"]
cell_cols  = [:firebrick, :darkorange, :forestgreen, :royalblue, :purple]

function global_ylims(rep_vals, xmax)
    vals = Float64[]
    for b6 in rep_vals
        sol = solutions_for_reps[b6]
        keep = findall(t -> t <= xmax, sol.t)
        for j in 1:5
            append!(vals, safe(sol[j,keep]))
        end
    end
    lo = max(minimum(vals), 1e-8)
    hi = maximum(vals) * 1.35
    return (lo, hi)
end

function make_trajectory_plot(rep_vals, xmax, filename; title_prefix="ODE trajectories")
    ylims_common = global_ylims(rep_vals, xmax)
    subplots = Any[]
    for (i, b6) in enumerate(rep_vals)
        sol = solutions_for_reps[b6]
        keep = findall(t -> t <= xmax, sol.t)
        t = sol.t[keep]

        p = plot(t, safe(sol[1,keep]),
                 yscale=:log10, ylims=ylims_common,
                 color=cell_cols[1], lw=2.5,
                 label=(i == 1 ? cell_names[1] : false),
                 xlabel="Time (days)", ylabel="Cells",
                 title="b6 = $(sci_label(b6))", legend=(i == 1 ? :topright : false),
                 margin=5Plots.mm)
        for j in 2:5
            plot!(p, t, safe(sol[j,keep]),
                  yscale=:log10, ylims=ylims_common,
                  color=cell_cols[j], lw=2.2,
                  label=(i == 1 ? cell_names[j] : false))
        end
        push!(subplots, p)
    end

    fig = plot(subplots..., layout=(length(rep_vals), 1), size=(1100, 220*length(rep_vals)))
    plot!(fig, plot_title=title_prefix)
    savefig(fig, joinpath(outdir, filename))
end


make_trajectory_plot(rep_global, 365.0, "b6_dense_trajectories_0_365_pretty.png",
                     title_prefix="Representative ODE trajectories (0–365 days)")
make_trajectory_plot(rep_global, 10000.0, "b6_dense_trajectories_0_10000_pretty.png",
                     title_prefix="Representative ODE trajectories (0–10,000 days)")

make_trajectory_plot(rep_thresh, 10000.0, "b6_threshold_trajectories_0_10000_pretty.png",
                     title_prefix="Threshold-region ODE trajectories (fine b6 grid, 0–10,000 days)")
println("Saved trajectory plots.")


# ODE summary
coarse_rows = rows[rows.grid_type .!= "fine", :]
fine_rows = rows[(rows.b6 .>= fine_lo) .& (rows.b6 .<= fine_hi), :]

function add_fine_band!(p)
    vspan!(p, [fine_lo, fine_hi], color=:gold, alpha=0.10, label=false)
end

p1 = plot(coarse_rows.b6, coarse_rows.T100, xscale=:log10, yscale=:log10,
          color=:lightgray, lw=1.6, label=false)
plot!(p1, coarse_rows.b6, coarse_rows.T365, xscale=:log10, yscale=:log10,
      color=:silver, lw=1.6, label=false)
plot!(p1, coarse_rows.b6, coarse_rows.T10000, xscale=:log10, yscale=:log10,
      color=:darkgray, lw=1.6, label=false)
plot!(p1, rows.b6, rows.T100, xscale=:log10, yscale=:log10,
      color=:dodgerblue, marker=:circle, lw=2.4, label="T100")
plot!(p1, rows.b6, rows.T365, marker=:circle, color=:orangered, lw=2.4, label="T365")
plot!(p1, rows.b6, rows.T10000, marker=:circle, color=:purple, lw=2.4, label="T10000")
add_fine_band!(p1)
hline!(p1, [1e3, 1e5], linestyle=:dash, color=:gray, label=false)
plot!(p1, xlabel="b6", ylabel="Tumor cells", title="Tumor output across the full b6 sweep")

p2 = plot(coarse_rows.b6, coarse_rows.t_to_1e5_plot, xscale=:log10,
          color=:lightgray, lw=1.6, label=false)
plot!(p2, coarse_rows.b6, coarse_rows.t_to_1e6_plot, xscale=:log10,
      color=:silver, lw=1.6, label=false)
plot!(p2, rows.b6, rows.t_to_1e5_plot, marker=:circle, color=:dodgerblue, lw=2.4, label="time to 1e5")
plot!(p2, rows.b6, rows.t_to_1e6_plot, marker=:circle, color=:orangered, lw=2.4, label="time to 1e6")
add_fine_band!(p2)
hline!(p2, [365.0, 10000.0], linestyle=:dash, color=:gray, label=false)
plot!(p2, xlabel="b6", ylabel="Days", title="Tumor threshold-crossing times")

p3 = plot(coarse_rows.b6, coarse_rows.auc_0_365, xscale=:log10, yscale=:log10,
          color=:lightgray, lw=1.6, label=false)
plot!(p3, coarse_rows.b6, coarse_rows.auc_0_1000, xscale=:log10, yscale=:log10,
      color=:silver, lw=1.6, label=false)
plot!(p3, rows.b6, rows.auc_0_365, marker=:circle, color=:dodgerblue, lw=2.4, label="AUC 0–365")
plot!(p3, rows.b6, rows.auc_0_1000, marker=:circle, color=:orangered, lw=2.4, label="AUC 0–1000")
add_fine_band!(p3)
plot!(p3, xlabel="b6", ylabel="Tumor AUC", title="Early tumor burden")

p4 = plot(rows.b6, rows.B10000, xscale=:log10, yscale=:log10,
          marker=:circle, color=:purple, lw=2.4, label="B")
plot!(p4, rows.b6, rows.NK10000, marker=:circle, color=:forestgreen, lw=2.4, label="NK")
plot!(p4, rows.b6, rows.CTL10000, marker=:circle, color=:royalblue, lw=2.4, label="CTL")
add_fine_band!(p4)
plot!(p4, xlabel="b6", ylabel="Cells at day 10,000", title="Immune compartments at long time")

summary_fig = plot(p1, p2, p3, p4, layout=(2,2), size=(1300, 950), margin=6Plots.mm)
savefig(summary_fig, joinpath(outdir, "b6_dense_ode_summary_pretty.png"))

#  zoomed threshold plot
q1 = plot(fine_rows.b6, fine_rows.T365, xscale=:log10, yscale=:log10,
          marker=:circle, color=:orangered, lw=2.6, label="T365")
plot!(q1, fine_rows.b6, fine_rows.T10000, marker=:circle, color=:purple, lw=2.6, label="T10000")
hline!(q1, [1e3, 1e5], linestyle=:dash, color=:gray, label=false)
plot!(q1, xlabel="b6", ylabel="Tumor cells", title="Zoomed threshold plot: T365 vs T10000")

q2 = plot(fine_rows.b6, fine_rows.t_to_1e5_plot, xscale=:log10,
          marker=:circle, color=:dodgerblue, lw=2.6, label="time to 1e5")
plot!(q2, fine_rows.b6, fine_rows.t_to_1e6_plot,
      marker=:circle, color=:orangered, lw=2.6, label="time to 1e6")
hline!(q2, [365.0, 10000.0], linestyle=:dash, color=:gray, label=false)
plot!(q2, xlabel="b6", ylabel="Days", title="Zoomed threshold plot: crossing times")

q3 = plot(fine_rows.b6, fine_rows.auc_0_365, xscale=:log10, yscale=:log10,
          marker=:circle, color=:dodgerblue, lw=2.6, label="AUC 0–365")
plot!(q3, fine_rows.b6, fine_rows.auc_0_1000,
      marker=:circle, color=:orangered, lw=2.6, label="AUC 0–1000")
plot!(q3, xlabel="b6", ylabel="Tumor AUC", title="Zoomed threshold plot: early burden")

# add a  regime-style panel using log10
regime_vals = log10.(fine_rows.T10000 .+ 1.0)
q4 = scatter(fine_rows.b6, regime_vals, xscale=:log10,
             markerstrokewidth=0, marker_z=regime_vals,
             color=:viridis, ms=9, label=false)
plot!(q4, fine_rows.b6, regime_vals, color=:black, alpha=0.35, lw=1.3, label=false)
plot!(q4, xlabel="b6", ylabel="log10(T10000 + 1)", title="Zoomed threshold plot: long-time state")

zoom_fig = plot(q1, q2, q3, q4, layout=(2,2), size=(1300, 950), margin=6Plots.mm)
savefig(zoom_fig, joinpath(outdir, "b6_dense_threshold_zoom_pretty.png"))
println("Saved summary + zoom plots.")

# console summary

controlled = rows[rows.final_state .== "tumor_controlled", :]
if nrow(controlled) > 0
    println("Smallest b6 classified as tumor-controlled by T(10000)<1e3: ", minimum(controlled.b6))
else
    println("No b6 value classified as tumor-controlled by T(10000)<1e3")
end

println("Done. Outputs are in: ", outdir)
