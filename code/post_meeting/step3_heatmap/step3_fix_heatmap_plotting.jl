using Plots
using DataFrames
using CSV
using Statistics


# Step 3 plotting fix only 

# Problem with the previous continuous heatmap:
#   The title said log10(T(10000)+1), but the contour overlay used raw T values,
#   so the colorbar got dominated/mislabeled by raw tumor counts.


# Fix:
#   1. Read the already-generated CSV.
#   2. Plot heatmaps in log-parameter coordinates: x = log10(b6), y = log10(b5).
#   3. Use logT10000 for the heatmap color with clims=(0,7).
#   4. Draw threshold boundaries as simple black/gray line overlays, not as a contour series
#      that changes the colorbar

outdir = "step3_fixed_plots"
mkpath(outdir)

csvfile = "step3_b5_b6_grid.csv"
if !isfile(csvfile)
    # fallback if running from a different working directory
    csvfile = joinpath(@__DIR__, "step3_b5_b6_grid.csv")
end

df = CSV.read(csvfile, DataFrame)

b5_vals = sort(unique(df.b5))
b6_vals = sort(unique(df.b6))
logb5 = log10.(b5_vals)
logb6 = log10.(b6_vals)

function build_matrix(df, colname, b5_values, b6_values)
    Z = fill(NaN, length(b5_values), length(b6_values))
    for (i, b5) in enumerate(b5_values)
        for (j, b6) in enumerate(b6_values)
            sub = df[(df.b5 .== b5) .& (df.b6 .== b6), :]
            if nrow(sub) == 1
                Z[i,j] = sub[1, colname]
            end
        end
    end
    return Z
end

Zlog = build_matrix(df, :logT10000, b5_vals, b6_vals)
Zphase = build_matrix(df, :phase_code, b5_vals, b6_vals)

# Boundary helper: smallest b6 where T(10000) < cutoff for each b5.
function boundary_for_cutoff(df, cutoff)
    xs = Float64[]
    ys = Float64[]
    for b5 in b5_vals
        sub = sort(df[df.b5 .== b5, :], :b6)
        idx = findfirst(sub.T10000 .< cutoff)
        if !isnothing(idx)
            push!(xs, log10(sub.b6[idx]))
            push!(ys, log10(b5))
        end
    end
    return xs, ys
end

bnd1_x, bnd1_y = boundary_for_cutoff(df, 1e3)
bnd2_x, bnd2_y = boundary_for_cutoff(df, 1e5)

# ------------------------------------------------------------
# Fixed continuous heatmap
# ------------------------------------------------------------
p1 = heatmap(
    logb6, logb5, Zlog,
    color=cgrad(:Blues),
    clims=(0,7),
    xlabel="log10(b6)",
    ylabel="log10(b5)",
    colorbar_title="log10(T(10000)+1)",
    title="Step 3: b5 × b6 long-time tumor burden",
    size=(900,700),
    right_margin=9Plots.mm,
    left_margin=8Plots.mm,
    bottom_margin=6Plots.mm,
    top_margin=6Plots.mm,
)
plot!(p1, bnd1_x, bnd1_y, color=:black, lw=3, label="T(10000)<1e3 boundary")
plot!(p1, bnd2_x, bnd2_y, color=:gray40, lw=2.5, linestyle=:dash, label="T(10000)<1e5 boundary")
savefig(p1, joinpath(outdir, "step3_fixed_log_heatmap.png"))

# ------------------------------------------------------------
# Fixed discrete phase map
# phase_code: 0=controlled, 1=intermediate, 2=tumor-dominant
# ------------------------------------------------------------
phase_colors = cgrad([:white, :lightskyblue, :navy])
p2 = heatmap(
    logb6, logb5, Zphase,
    color=phase_colors,
    clims=(0,2),
    xlabel="log10(b6)",
    ylabel="log10(b5)",
    colorbar_ticks=([0,1,2], ["controlled", "intermediate", "dominant"]),
    title="Step 3: discrete long-time regimes",
    size=(900,700),
    right_margin=12Plots.mm,
    left_margin=8Plots.mm,
    bottom_margin=6Plots.mm,
    top_margin=6Plots.mm,
)
plot!(p2, bnd1_x, bnd1_y, color=:black, lw=3, label="T(10000)<1e3 boundary")
savefig(p2, joinpath(outdir, "step3_fixed_phase_map.png"))

# ------------------------------------------------------------
# Boundary plot: this is the cleanest figure for the story.
# ------------------------------------------------------------
# Convert boundary points back to b-scale for a more readable log-log line plot.
boundary_df = DataFrame(
    log10_b5 = bnd1_y,
    log10_b6_control = bnd1_x,
    b5 = 10 .^ bnd1_y,
    b6_control_boundary = 10 .^ bnd1_x,
)
CSV.write(joinpath(outdir, "step3_fixed_control_boundary.csv"), boundary_df)

p3 = plot(
    boundary_df.b5, boundary_df.b6_control_boundary,
    marker=:circle,
    lw=3,
    xscale=:log10,
    yscale=:log10,
    xlabel="b5",
    ylabel="b6 needed for T(10000)<1e3",
    title="Step 3: stronger b5 requires larger b6 for control",
    legend=false,
    size=(850,550),
    left_margin=8Plots.mm,
    bottom_margin=6Plots.mm,
)
savefig(p3, joinpath(outdir, "step3_fixed_boundary_plot.png"))

println("Fixed Step 3 plots saved in: ", outdir)
println("Main figures to show:")
println("  1) step3_fixed_log_heatmap.png")
println("  2) step3_fixed_phase_map.png")
println("  3) step3_fixed_boundary_plot.png")
