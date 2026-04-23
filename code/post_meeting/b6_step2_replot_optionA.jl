using CSV
using DataFrames
using DifferentialEquations
using Plots
using Printf

# =========================================================
# Step 2 modified v - replot existing initial-condition robustness results
# Goal: make the plots clearer without rerunning the full 1000-run sweep.


ENV["GKSwstype"] = "100"
gr()

script_dir = @__DIR__
outdir = joinpath(script_dir, "b6_ic_robustness_replots")
mkpath(outdir)

function find_input_file(fname)
    candidate_paths = unique([
        abspath(joinpath(script_dir, fname)),
        abspath(joinpath(script_dir, "b6_ic_robustness_outputs", fname)),
        abspath(joinpath(script_dir, "..", "d1", "b6_ic_robustness_outputs", fname)),
        abspath(fname),
        abspath(joinpath("b6_ic_robustness_outputs", fname)),
        joinpath("/mnt/data", fname)
    ])

    for path in candidate_paths
        if isfile(path)
            println("Using input file: ", path)
            return path
        end
    end

    error("Cannot find input file: $(fname)\nSearched:\n  " * join(candidate_paths, "\n  "))
end

all_runs = CSV.read(find_input_file("b6_ic_robustness_all_runs.csv"), DataFrame)
thresholds = CSV.read(find_input_file("b6_ic_robustness_thresholds.csv"), DataFrame)

changed_comps = ["Tumor", "MDSC", "NK", "CTL", "B"]
mults = [0.2, 0.5, 1.0, 2.0, 5.0]
comp_colors = Dict("Tumor"=>:red, "MDSC"=>:orange, "NK"=>:green, "CTL"=>:blue, "B"=>:purple)
comp_markers = Dict("Tumor"=>:circle, "MDSC"=>:rect, "NK"=>:utriangle, "CTL"=>:diamond, "B"=>:star5)
traj_b6 = 8.0e-7

safe(v; eps=1e-8) = max.(collect(v), eps)

function padded_limits(values; pad_frac=0.08, min_pad=1e-8)
    finite_vals = [v for v in values if isfinite(v)]
    isempty(finite_vals) && return (0.0, 1.0)

    lo, hi = extrema(finite_vals)
    pad = max((hi - lo) * pad_frac, min_pad)
    if hi == lo
        pad = max(abs(lo) * pad_frac, min_pad)
    end
    return (lo - pad, hi + pad)
end

format_sci(x) = @sprintf("%.2e", x)

function spread_label(values)
    vals = [v for v in values if isfinite(v) && v > 0.0]
    isempty(vals) && return "NA"

    spread = maximum(vals) / minimum(vals)
    if spread >= 100.0
        return @sprintf("%.0fx", spread)
    elseif spread >= 10.0
        return @sprintf("%.1fx", spread)
    else
        return @sprintf("%.2fx", spread)
    end
end

function threshold_shift_info(col::Symbol)
    best = (comp="none", mult=1.0, rel=1.0, pct=0.0)
    max_abs_shift = -Inf

    for comp in changed_comps
        sub = thresholds[thresholds.changed_comp .== comp, :]
        sort!(sub, :multiplier)
        base = sub[sub.multiplier .== 1.0, col][1]

        for row in eachrow(sub)
            rel = row[col] / base
            pct = 100.0 * (rel - 1.0)
            if abs(pct) > max_abs_shift
                max_abs_shift = abs(pct)
                best = (comp=comp, mult=row.multiplier, rel=rel, pct=pct)
            end
        end
    end

    return best
end

function threshold_note(info; quiet_cutoff=0.5)
    if abs(info.pct) < quiet_cutoff
        return "No visible threshold shift across IC perturbations"
    end
    return @sprintf("Largest shift: %s(0) at %.1fx (%+.1f%%)",
                    info.comp, info.mult, info.pct)
end

function add_panel_note!(p, note; x=0.24, y=0.95)
    xl = Plots.xlims(p)
    yl = Plots.ylims(p)
    xpos = xl[1] + x * (xl[2] - xl[1])
    ypos = yl[1] + y * (yl[2] - yl[1])
    annotate!(p, xpos, ypos, text(note, 9, :left, :black))
end

# Plot 1A: relative threshold plot

function make_relative_threshold_plot()
    p1 = plot(title="Relative threshold for T(10000) < 1e5",
              xlabel="Initial-condition multiplier",
              ylabel="Threshold / baseline threshold",
              xscale=:log10,
              legend=:outerright,
              size=(1200,450),
              gridalpha=0.25)
    p2 = plot(title="Relative threshold for T(10000) < 1e3",
              xlabel="Initial-condition multiplier",
              ylabel="Threshold / baseline threshold",
              xscale=:log10,
              legend=false,
              gridalpha=0.25)

    hspan!(p1, [0.98, 1.02], color=:gray90, alpha=0.8, label=false)
    hspan!(p2, [0.98, 1.02], color=:gray90, alpha=0.8, label=false)
    hline!(p1, [1.0], color=:gray, linestyle=:dash, lw=2, label="baseline")
    hline!(p2, [1.0], color=:gray, linestyle=:dash, lw=2, label="baseline")
    rel1_all = Float64[]
    rel3_all = Float64[]

    for comp in changed_comps
        sub = thresholds[thresholds.changed_comp .== comp, :]
        sort!(sub, :multiplier)

        base1 = sub[sub.multiplier .== 1.0, :threshold_T1e5][1]
        base3 = sub[sub.multiplier .== 1.0, :threshold_T1e3][1]

        rel1 = sub.threshold_T1e5 ./ base1
        rel3 = sub.threshold_T1e3 ./ base3
        append!(rel1_all, rel1)
        append!(rel3_all, rel3)

        plot!(p1, sub.multiplier, rel1,
              marker=comp_markers[comp], lw=2.5, ms=6,
              color=comp_colors[comp], label=comp)
        plot!(p2, sub.multiplier, rel3,
              marker=comp_markers[comp], lw=2.5, ms=6,
              color=comp_colors[comp], label=comp)
    end

    y1 = padded_limits(rel1_all; pad_frac=0.18, min_pad=0.02)
    y3 = padded_limits(rel3_all; pad_frac=0.18, min_pad=0.02)

    plot!(p1, xticks=(mults, string.(mults)), ylim=(max(0.85, y1[1]), y1[2]))
    plot!(p2, xticks=(mults, string.(mults)), ylim=(max(0.90, y3[1]), y3[2]))

    add_panel_note!(p1, threshold_note(threshold_shift_info(:threshold_T1e5)); x=0.10, y=0.93)
    add_panel_note!(p2, threshold_note(threshold_shift_info(:threshold_T1e3)); x=0.10, y=0.93)

    fig = plot(p1, p2, layout=(1,2), size=(1250,460), margin=5Plots.mm)
    savefig(fig, joinpath(outdir, "ic_thresholds_relative_better.png"))
end
make_relative_threshold_plot()

# Plot 1B: raw thresholds, useful as backup
function make_raw_threshold_plot()
    p1 = plot(title="Raw threshold: T(10000) < 1e5",
              xlabel="Initial-condition multiplier", ylabel="b6 threshold",
              xscale=:log10, yscale=:identity,
              yformatter=format_sci, legend=:outerright, gridalpha=0.25)
    p2 = plot(title="Raw threshold: T(10000) < 1e3",
              xlabel="Initial-condition multiplier", ylabel="b6 threshold",
              xscale=:log10, yscale=:identity,
              yformatter=format_sci, legend=false, gridalpha=0.25)

    raw1_all = Float64[]
    raw3_all = Float64[]
    base_row = thresholds[(thresholds.changed_comp .== "Tumor") .& (thresholds.multiplier .== 1.0), :]
    base1 = base_row.threshold_T1e5[1]
    base3 = base_row.threshold_T1e3[1]

    hline!(p1, [base1], color=:gray50, linestyle=:dash, lw=2, label="baseline")
    hline!(p2, [base3], color=:gray50, linestyle=:dash, lw=2, label="baseline")

    for comp in changed_comps
        sub = thresholds[thresholds.changed_comp .== comp, :]
        sort!(sub, :multiplier)
        append!(raw1_all, sub.threshold_T1e5)
        append!(raw3_all, sub.threshold_T1e3)
        plot!(p1, sub.multiplier, sub.threshold_T1e5,
              marker=comp_markers[comp], lw=2.5, ms=6,
              color=comp_colors[comp], label=comp)
        plot!(p2, sub.multiplier, sub.threshold_T1e3,
              marker=comp_markers[comp], lw=2.5, ms=6,
              color=comp_colors[comp], label=comp)
    end

    plot!(p1,
          xticks=(mults, string.(mults)),
          ylim=padded_limits(raw1_all; pad_frac=0.20, min_pad=2e-8))
    plot!(p2,
          xticks=(mults, string.(mults)),
          ylim=padded_limits(raw3_all; pad_frac=0.12, min_pad=2e-8))

    add_panel_note!(p1, threshold_note(threshold_shift_info(:threshold_T1e5)); x=0.08, y=0.92)
    add_panel_note!(p2, threshold_note(threshold_shift_info(:threshold_T1e3)); x=0.08, y=0.92)

    fig = plot(p1, p2, layout=(1,2), size=(1250,460), margin=5Plots.mm)
    savefig(fig, joinpath(outdir, "ic_thresholds_raw_better.png"))
end
make_raw_threshold_plot()

# Plot 2: discrete phase heatmap
# This is easier to interpret than a continuous heatmap because the question is regime-level
# Categories:
#   controlled = T(10000) < 1e3
#   intermediate = 1e3 <= T(10000) < 1e5
#   tumor_dominant = T(10000) >= 1e5
function phase_number(state)
    if state == "controlled"
        return 1.0
    elseif state == "intermediate"
        return 2.0
    else
        return 3.0
    end
end

function make_phase_heatmap()
    b6_grid = sort(unique(all_runs.b6))
    logb6 = log10.(b6_grid)
    hm_plots = Any[]

    for comp in changed_comps
        Z = fill(NaN, length(mults), length(b6_grid))
        for (i, m) in enumerate(mults)
            for (j, b6) in enumerate(b6_grid)
                row = all_runs[(all_runs.changed_comp .== comp) .&
                               (all_runs.multiplier .== m) .&
                               (abs.(all_runs.b6 .- b6) .< 1e-20), :]
                if nrow(row) > 0
                    Z[i,j] = phase_number(row.final_state[1])
                end
            end
        end

        p = heatmap(logb6, mults, Z,
                    title="Change $(comp)(0)",
                    xlabel="log10(b6)", ylabel="IC multiplier",
                    color=cgrad([:darkseagreen2, :khaki1, :firebrick2], categorical=true),
                    clims=(1,3),
                    colorbar=false,
                    yticks=(mults, string.(mults)))

        sub = thresholds[thresholds.changed_comp .== comp, :]
        sort!(sub, :multiplier)
        plot!(p, log10.(sub.threshold_T1e5), sub.multiplier,
              color=:white, lw=2.2, marker=:circle, ms=4, label=false)
        plot!(p, log10.(sub.threshold_T1e3), sub.multiplier,
              color=:black, lw=2.0, linestyle=:dash, marker=:diamond, ms=4, label=false)
        push!(hm_plots, p)
    end

    legend_panel = plot(framestyle=:none, grid=false, axis=nothing,
                        legend=:best, title="Regime / threshold key")
    scatter!(legend_panel, [NaN], [NaN], marker=:square, ms=10,
             color=:darkseagreen2, label="controlled")
    scatter!(legend_panel, [NaN], [NaN], marker=:square, ms=10,
             color=:khaki1, label="intermediate")
    scatter!(legend_panel, [NaN], [NaN], marker=:square, ms=10,
             color=:firebrick2, label="tumor-dominant")
    plot!(legend_panel, [NaN], [NaN], color=:gray20, lw=2.2,
          marker=:circle, ms=4, label="T(10000) < 1e5")
    plot!(legend_panel, [NaN], [NaN], color=:gray20, lw=2.0,
          linestyle=:dash, marker=:diamond, ms=4, label="T(10000) < 1e3")

    fig = plot(hm_plots..., legend_panel, layout=(3,2), size=(1250,1100), margin=4Plots.mm)
    savefig(fig, joinpath(outdir, "ic_discrete_phase_heatmap.png"))
end
make_phase_heatmap()

# Plot 3: cleaner tumor trajectory replots near the threshold
# These rerun only 25 small ODE trajectories because the CSV stores summaries, not full paths
function base_p(; b6=traj_b6)
    return (
        a1=1.0e-1, h=1.0e7, b1=3.5e-6, b2=1.1e-7, z1=0.0,
        a2=1.0e2, a3=1.0e8, z2=2.0e-1,
        a4=1.4e4, a5=2.5e-2, b3=4.0e-5, z3=4.12e-2,
        a6=1.1e-7, a7=1.0e-1, b4=1.0e-4, z4=2.0e-2,
        g1=1.0e10, g2=2.02e7, g3=2.02e7,
        a8=1.0e4, a9=5.0e-2, b5=1.0e-4, z5=2.0e-2,
        g4=2.02e7, g5=1.0e3,
        a10=2.0, a11=1.0e-1, b6=b6
    )
end

function baseline_state(p; tumor_cells=2.0)
    xM0 = p.a2 / p.z2
    xN0 = p.z2 * p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC0 = 0.0
    xB0 = p.a8 / (p.z5 + p.b5 * xM0)
    return [tumor_cells, xM0, xN0, xC0, xB0]
end

function start_state_with_multiplier(p, changed_comp::String, mult::Float64; tumor_cells=2.0)
    u0 = baseline_state(p; tumor_cells=tumor_cells)
    comp_index = Dict("Tumor"=>1, "MDSC"=>2, "NK"=>3, "CTL"=>4, "B"=>5)
    idx = comp_index[changed_comp]

    if changed_comp == "CTL" && u0[idx] == 0.0
        ctl_ref = 1654.4182572290226
        u0[idx] = mult * ctl_ref
    else
        u0[idx] *= mult
    end
    return u0
end

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

function solve_case(b6, changed_comp, mult; tmax=10000.0, saveat=10.0)
    p = base_p(b6=b6)
    u0 = start_state_with_multiplier(p, changed_comp, mult)
    prob = ODEProblem(model!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5();
                saveat=saveat,
                reltol=1e-7,
                abstol=1e-9,
                callback=cb_bad,
                isoutofdomain=(u,p,t)->any(x->x < -1e-6, u))
    return sol
end

function line_style_for_multiplier(mult)
    if mult < 1.0
        return :dash
    elseif mult > 1.0
        return :solid
    end
    return :solid
end

function line_width_for_multiplier(mult)
    return mult == 1.0 ? 3.2 : 2.4
end

function line_alpha_for_multiplier(mult)
    return mult == 1.0 ? 1.0 : 0.9
end

function make_tumor_trajectory_replot(tmax, saveat, fname; log_time=false, start_time=0.0, title_text="")
    colors = [:black, :gray40, :dodgerblue3, :darkorange2, :crimson]
    subplots = Any[]

    for comp in changed_comps
        p = plot(xlabel="Time (days)", ylabel="Tumor cells",
                 title="Vary $(comp)(0)",
                 yscale=:log10,
                 xscale=log_time ? :log10 : :identity,
                 legend=false,
                 gridalpha=0.25)
        end_values = Float64[]

        for (i, mult) in enumerate(mults)
            sol = solve_case(traj_b6, comp, mult; tmax=tmax, saveat=saveat)
            times = collect(sol.t)
            values = safe(sol[1,:])

            if log_time
                keep = times .>= start_time
                times = times[keep]
                values = values[keep]
            end

            plot!(p, times, values,
                  lw=line_width_for_multiplier(mult),
                  ls=line_style_for_multiplier(mult),
                  alpha=line_alpha_for_multiplier(mult),
                  color=colors[i], label="$(mult)x")
            push!(end_values, last(values))
        end

        if log_time
            vline!(p, [365.0], color=:gray45, linestyle=:dot, lw=1.7, label=false)
            plot!(p, xlim=(start_time, tmax), ylim=(1.0, 2.0e7))
            annotate!(p, start_time * 1.15, 9.0e6,
                      text("spread @ $(Int(round(tmax)))d: $(spread_label(end_values))", 9, :left, :black))
        else
            plot!(p, xlim=(0.0, tmax), ylim=(1.0, 2.0e7))
            annotate!(p, 0.04 * tmax, 9.0e6,
                      text("spread @ $(Int(round(tmax)))d: $(spread_label(end_values))", 9, :left, :black))
        end
        push!(subplots, p)
    end

    legend_panel = plot(framestyle=:none, grid=false, axis=nothing,
                        legend=:best, title="IC multiplier")
    for (i, mult) in enumerate(mults)
        plot!(legend_panel, [NaN], [NaN], color=colors[i], lw=2.4, label="$(mult)x")
    end

    fig = plot(subplots..., legend_panel,
               layout=(3,2),
               size=(1280,1100),
               margin=4Plots.mm,
               plot_title=title_text)
    savefig(fig, joinpath(outdir, fname))
end

function make_nk_focus_trajectory()
    colors = [:black, :gray40, :dodgerblue3, :darkorange2, :crimson]
    p1 = plot(title="NK(0) perturbation: early-time divergence",
              xlabel="Time (days)", ylabel="Tumor cells",
              yscale=:log10, legend=:outerright, gridalpha=0.25)
    p2 = plot(title="NK(0) perturbation: long-time fate",
              xlabel="Time (days, log scale)", ylabel="Tumor cells",
              xscale=:log10, yscale=:log10, legend=false, gridalpha=0.25)

    long_end_values = Float64[]
    for (i, mult) in enumerate(mults)
        sol = solve_case(traj_b6, "NK", mult; tmax=10000.0, saveat=1.0)
        plot!(p1, sol.t[sol.t .<= 365.0], safe(sol[1, sol.t .<= 365.0]),
              lw=line_width_for_multiplier(mult),
              ls=line_style_for_multiplier(mult),
              alpha=line_alpha_for_multiplier(mult),
              color=colors[i], label="$(mult)x")

        keep = sol.t .>= 1.0
        plot!(p2, sol.t[keep], safe(sol[1, keep]),
              lw=line_width_for_multiplier(mult),
              ls=line_style_for_multiplier(mult),
              alpha=line_alpha_for_multiplier(mult),
              color=colors[i], label=false)
        push!(long_end_values, safe(sol[1, keep])[end])
    end

    plot!(p1, xlim=(0.0, 365.0), ylim=(1.0, 2.0e7))
    plot!(p2, xlim=(1.0, 10000.0), ylim=(1.0, 2.0e7))
    vline!(p2, [365.0], color=:gray45, linestyle=:dot, lw=1.7, label=false)

    annotate!(p1, 15.0, 9.0e6,
              text("Low NK(0) delays control and can re-enter escape", 10, :left, :black))
    annotate!(p2, 1.3, 9.0e6,
              text("spread @ 10000d: $(spread_label(long_end_values))", 10, :left, :black))

    fig = plot(p1, p2, layout=(1,2), size=(1280,500), margin=5Plots.mm,
               plot_title="Step 2 highlight: NK(0) is the only IC perturbation with a strong effect near threshold")
    savefig(fig, joinpath(outdir, "ic_nk_focus_trajectories_b6_8em7.png"))
end

make_tumor_trajectory_replot(365.0, 1.0,
                             "ic_tumor_trajectories_0_365_b6_8em7.png";
                             title_text="Near-threshold trajectories at b6 = 8e-7 (0-365 days)")
make_tumor_trajectory_replot(10000.0, 10.0,
                             "ic_tumor_trajectories_0_10000_b6_8em7.png";
                             log_time=true,
                             start_time=1.0,
                             title_text="Near-threshold trajectories at b6 = 8e-7 (1-10000 days, log-time)")
make_nk_focus_trajectory()

println("Done. Replots saved to: ", outdir)
