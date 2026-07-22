using CSV
using DataFrames
using Plots
using Statistics

# ------------------------------------------------------------
# Basin-of-attraction postprocessing / clean figures
# ------------------------------------------------------------
# This script reads the basin-of-attraction outputs and recreates the
# cleaner summary figures used for interpretation.
#
# It does NOT rerun the ODE simulations.
#
# Expected input:
#   basin_of_attraction_outputs/basin_all_runs.csv
#   basin_of_attraction_outputs/basin_fraction_summary.csv
#
# Output folder:
#   basin_interpretation_figures/
#
# Figure-title rule:
#   Use descriptive result-based titles only.
#   No "Phase 1", "Step 4", etc.
# ------------------------------------------------------------

outdir = "basin_interpretation_figures"
mkpath(outdir)

function first_existing(paths)
    for p in paths
        if isfile(p)
            return p
        end
    end
    error("Could not find expected input file in any of these paths:\n" * join(paths, "\n"))
end

all_runs_file = first_existing([
    joinpath("basin_of_attraction_outputs", "basin_all_runs.csv"),
    "basin_all_runs.csv",
    joinpath("results", "post_meeting", "basin_of_attraction", "basin_all_runs.csv")
])

frac_file = first_existing([
    joinpath("basin_of_attraction_outputs", "basin_fraction_summary.csv"),
    "basin_fraction_summary.csv",
    joinpath("results", "post_meeting", "basin_of_attraction", "basin_fraction_summary.csv")
])

rows = CSV.read(all_runs_file, DataFrame)
frac = CSV.read(frac_file, DataFrame)

println("Using basin runs: ", all_runs_file)
println("Using basin fraction summary: ", frac_file)

# Regime color scheme:
# controlled = green
# intermediate = yellow
# tumor-dominant = red
phase_colors = cgrad([:forestgreen, :khaki, :firebrick], categorical=true)

function build_matrix(df, x_col, y_col, x_vals, y_vals, value_col)
    Z = fill(NaN, length(y_vals), length(x_vals))
    for (i, yv) in enumerate(y_vals)
        for (j, xv) in enumerate(x_vals)
            sub = df[(df[!, x_col] .== xv) .& (df[!, y_col] .== yv), :]
            if nrow(sub) == 1
                Z[i, j] = sub[1, value_col]
            end
        end
    end
    return Z
end

function label_map_name(map_name)
    if map_name == "Tumor0_NK0"
        return "Tumor(0) × NK(0)"
    elseif map_name == "Tumor0_B0"
        return "Tumor(0) × B(0)"
    elseif map_name == "NK0_B0"
        return "B(0) × NK(0)"
    else
        return map_name
    end
end

function axis_info(map_name)
    if map_name == "Tumor0_NK0"
        return (:tumor0, :nk0, "log10(Tumor(0))", "log10(NK(0))")
    elseif map_name == "Tumor0_B0"
        return (:tumor0, :b0, "log10(Tumor(0))", "log10(B(0))")
    elseif map_name == "NK0_B0"
        return (:b0, :nk0, "log10(B(0))", "log10(NK(0))")
    else
        error("Unknown map name: $(map_name)")
    end
end

function one_basin_panel(df, map_name, b6; title="")
    x_col, y_col, xlab, ylab = axis_info(map_name)
    sub = df[(df.map_name .== map_name) .& (df.b6 .== b6), :]
    x_vals = sort(unique(sub[!, x_col]))
    y_vals = sort(unique(sub[!, y_col]))
    Z = build_matrix(sub, x_col, y_col, x_vals, y_vals, :regime_code)

    p = heatmap(
        log10.(x_vals), log10.(y_vals), Z,
        color=phase_colors,
        clims=(0,2),
        xlabel=xlab,
        ylabel=ylab,
        title=title,
        colorbar=false,
        size=(550, 450),
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
        top_margin=6Plots.mm
    )
    return p
end

# ------------------------------------------------------------
# Figure A: Tumor(0) × NK(0) maps across b6 values
# ------------------------------------------------------------
b6_values = sort(unique(rows.b6))
map_name = "Tumor0_NK0"

panels = []
for b6 in b6_values
    push!(panels, one_basin_panel(rows, map_name, b6; title="b6 = $(round(b6, sigdigits=3))"))
end

p_tn = plot(
    panels...,
    layout=(1, length(panels)),
    size=(520*length(panels), 460),
    plot_title="Higher b6 expands the tumor-control basin in the Tumor(0) × NK(0) plane"
)
savefig(p_tn, joinpath(outdir, "basin_tumor_nk_across_b6.png"))

# ------------------------------------------------------------
# Figure B: compare all three maps at b6 = 1.5e-6 (or nearest)
# ------------------------------------------------------------
target_b6 = 1.5e-6
b6_nearest = b6_values[argmin(abs.(b6_values .- target_b6))]

map_names = ["Tumor0_NK0", "Tumor0_B0", "NK0_B0"]
panels = []
for mn in map_names
    push!(panels, one_basin_panel(rows, mn, b6_nearest; title=label_map_name(mn)))
end

p_three = plot(
    panels...,
    layout=(1,3),
    size=(1500, 470),
    plot_title="At the same b6, initial immune context changes long-time tumor fate"
)
savefig(p_three, joinpath(outdir, "basin_three_maps_b6_1p5e-6.png"))

# ------------------------------------------------------------
# Figure C: controlled basin fraction vs b6
# ------------------------------------------------------------
p_frac = plot(
    xlabel="b6",
    ylabel="fraction of IC grid ending controlled",
    xscale=:log10,
    title="The controlled basin expands as b6 increases",
    size=(850, 550),
    legend=:topleft,
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)

for mn in unique(frac.map_name)
    sub = frac[frac.map_name .== mn, :]
    sub = sort(sub, :b6)
    plot!(p_frac, sub.b6, sub.frac_controlled,
          marker=:circle, lw=2.8, label=label_map_name(mn))
end

savefig(p_frac, joinpath(outdir, "controlled_basin_fraction_clean.png"))

# ------------------------------------------------------------
# Figure D: stacked final-regime fraction summary
# ------------------------------------------------------------
# This is useful as a supplement. It shows all regimes, not only controlled.
labels = String[]
controlled = Float64[]
intermediate = Float64[]
dominant = Float64[]

for mn in map_names
    for b6 in b6_values
        sub = frac[(frac.map_name .== mn) .& (frac.b6 .== b6), :]
        if nrow(sub) == 1
            push!(labels, replace(label_map_name(mn), " × " => "\n× ") * "\n" * string(round(b6, sigdigits=3)))
            push!(controlled, sub.frac_controlled[1])
            push!(intermediate, sub.frac_intermediate[1])
            push!(dominant, sub.frac_dominant[1])
        end
    end
end

x = 1:length(labels)
p_stack = bar(
    x, controlled,
    label="controlled",
    color=:forestgreen,
    ylabel="fraction of grid",
    xticks=(x, labels),
    xrotation=45,
    title="Final-regime fractions across basin maps",
    size=(1100, 600),
    legend=:topright,
    bottom_margin=18Plots.mm,
    left_margin=8Plots.mm,
    top_margin=8Plots.mm
)
bar!(p_stack, x, intermediate, bottom=controlled, label="intermediate", color=:khaki)
bar!(p_stack, x, dominant, bottom=controlled .+ intermediate, label="tumor-dominant", color=:firebrick)

savefig(p_stack, joinpath(outdir, "basin_regime_fraction_stacked.png"))

println("Clean basin interpretation figures saved in: ", outdir)
println("Main figures:")
println("  basin_tumor_nk_across_b6.png")
println("  basin_three_maps_b6_1p5e-6.png")
println("  controlled_basin_fraction_clean.png")
println("Supplement:")
println("  basin_regime_fraction_stacked.png")
