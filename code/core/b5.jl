using DifferentialEquations
using Plots
using Colors
using Random
using Statistics
using DataFrames
using CSV

# B-cell extension of the Kreger model


logo = Colors.JULIA_LOGO_COLORS
cell_colors = [logo.red, colorant"goldenrod2", logo.green, logo.blue, colorant"purple"]

# baseline parameters

function base_params(; b5 = 1.0e-4, a10 = 2.0, a11 = 0.0, b6 = 0.0, ctl_mode = :tumor_gated)
    return (
        # original CIR parameters
        a1 = 1.0e-1,
        h= 1.0e7,
        b1 = 3.5e-6,
        b2 = 1.1e-7,
        z1 = 0.0,

        a2 = 1.0e2,
        a3 = 1.0e8,
        z2 = 2.0e-1,

        a4 = 1.4e4,
        a5 = 2.5e-2,
        b3= 4.0e-5,
        z3= 4.12e-2,

        a6 = 1.1e-7,
        a7 = 1.0e-1,
        b4 = 1.0e-4,
        z4 = 2.0e-2,

        g1 = 1.0e10,
        g2 = 2.02e7,
        g3 = 2.02e7,

        #new B-cell parameters
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

# 
# initial states
function tumor_free_state(p;tumor_cells = 2.0)
    xM_star = p.a2 / p.z2
    xN_star = p.z2 * p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC_star = 0.0
    xB_star = p.a8 / (p.z5 + p.b5 * xM_star)
    return [tumor_cells, xM_star, xN_star, xC_star, xB_star]
end

#based on Kreger model
function day100_state(p)
    xM0 = 804.0710624094437
    xB0 = p.a8 / (p.z5 + p.b5 * xM0)
    return [8395.368084462818,
            804.0710624094437,
            197565.74910954377,
            1654.4182572290226,
            xB0]
end

# model equations
function model_B!(du, u, p, t)
    xT, xM, xN, xC, xB = u

    B_to_CTL = if p.ctl_mode == :direct
        p.a11 * (xB / (p.g5 + xB))
    else
        p.a11 * (xT^2 / (p.g3 + xT^2)) * (xB / (p.g5 + xB))
    end

    # tumor
    du[1] = p.a1 * xT * log(max(p.h / xT, 1.0)) - p.b1 * xT * xN - p.b2 * xT * xC - p.b6 * xT * xB - p.z1 * xT

    # MDSC
    du[2] = p.a2 + p.a3 * xT / (p.g1 + xT) - p.z2 * xM

    # NK
    du[3] = p.a4 + p.a5 * (xT^2) / (p.g2 + xT^2) * (1.0 + p.a10 * xB / (p.g5 + xB)) - p.b3 * xM * xN - p.z3 * xN

    # CTL
    du[4] = p.a6 * xT * xN + p.a7 * (xT^2) / (p.g3 + xT^2) + B_to_CTL - p.b4 * xM * xC - p.z4 * xC

    # effector B
    du[5] = p.a8 + p.a9 * (xT^2) / (p.g4 + xT^2) - p.b5 * xM * xB - p.z5 * xB
end

function sigma_B!(du, u, p, t)
    for i in eachindex(u)
        du[i] = u[i]
    end
end

function negative_condition(u, t, integrator)
    any(x -> x < 0.0, u)
end
stop_if_negative = DiscreteCallback(negative_condition, terminate!)

#solve helpers
function solve_ode_b5(b5; a10 = 2.0, a11 = 0.0, b6 = 0.0,
                      ctl_mode = :tumor_gated, init_mode = :tumor_free,
                      tmax = 365.0, saveat = 1.0)
    p = base_params(b5 = b5, a10 = a10, a11 = a11, b6 = b6, ctl_mode = ctl_mode)
    u0 = init_mode == :day100 ? day100_state(p) : tumor_free_state(p; tumor_cells = 2.0)
    prob = ODEProblem(model_B!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5(); saveat = saveat,
                reltol = 1e-6, abstol = 1e-9,
                callback = stop_if_negative,
                isoutofdomain = (u,p,t)->any(x->x<0,u))
    return sol, p
end

#trajectory plot
function safe_log_vec(v; eps = 1e-6)
    return max.(collect(v), eps)
end

function plot_b5_trajectories(b5_values; a10 = 2.0, a11 = 0.0, b6 = 0.0,
                                    ctl_mode = :tumor_gated, init_mode = :tumor_free)
    subplots = Plots.Plot[]
    for b5 in b5_values
        sol, _ = solve_ode_b5(b5; a10 = a10, a11 = a11, b6 = b6,
                              ctl_mode = ctl_mode, init_mode = init_mode)

        p = plot(sol.t, safe_log_vec(sol[1,:]), lw = 3, yscale = :log10,
                 color = cell_colors[1], label = "Tumor", xlabel = "Time (days)",
                 ylabel = "Cells", title = "ODE trajectories, b5 = $(b5), a10 = $(a10)")
        plot!(p, sol.t, safe_log_vec(sol[2,:]), lw = 3, color = cell_colors[2], label = "MDSC")
        plot!(p, sol.t, safe_log_vec(sol[3,:]), lw = 3, color = cell_colors[3], label = "NK")
        plot!(p, sol.t, safe_log_vec(sol[4,:]), lw = 3, color = cell_colors[4], label = "CTL")
        plot!(p, sol.t, safe_log_vec(sol[5,:]), lw = 3, color = cell_colors[5], label = "B")
        push!(subplots, p)
    end
    return plot(subplots..., layout = (length(b5_values), 1), size = (1000, 260 * length(b5_values)))
end

#ODE sweep
function sweep_b5_ode(b5_values; a10 = 2.0, a11 = 0.0, b6 = 0.0,
                      ctl_mode = :tumor_gated, init_mode = :tumor_free,
                      outfile = "beta5_ode_sweep.csv")
    rows = DataFrame(b5 = Float64[], T365 = Float64[], M365 = Float64[], NK365 = Float64[], CTL365 = Float64[], B365 = Float64[])
    for b5 in b5_values
        sol, _ = solve_ode_b5(b5; a10 = a10, a11 = a11, b6 = b6,
                              ctl_mode = ctl_mode, init_mode = init_mode)
        push!(rows, (b5, sol[1,end], sol[2,end], sol[3,end], sol[4,end], sol[5,end]))
    end
    CSV.write(outfile, rows)
    return rows
end

function plot_b5_summary(df)
    p1 = plot(df.b5, df.T365, marker = :circle, lw = 3, xscale = :log10,
              xlabel = "b5", ylabel = "Tumor at day 365", title = "Tumor burden vs b5")
    p2 = plot(df.b5, df.B365, marker = :circle, lw = 3, xscale = :log10,
              xlabel = "b5", ylabel = "B at day 365", title = "B abundance vs b5")
    p3 = plot(df.b5, df.NK365, marker = :circle, lw = 3, xscale = :log10,
              xlabel = "b5", ylabel = "NK at day 365", title = "NK abundance vs b5")
    p4 = plot(df.b5, df.CTL365, marker = :circle, lw = 3, xscale = :log10,
              xlabel = "b5", ylabel = "CTL at day 365", title = "CTL abundance vs b5")
    return plot(p1, p2, p3, p4, layout = (2,2), size = (1000, 800))
end

#SDE
function stochastic_b5_summary(b5_values; a10 = 2.0, a11 = 0.0, b6 = 0.0,
                               ctl_mode = :tumor_gated, nreps = 100, seed = 14,
                               outfile = "beta5_sde_summary.csv")
    Random.seed!(seed)
    out = DataFrame(b5 = Float64[], success_prob = Float64[], mean_success_tumor = Float64[], mean_time_to_extinction = Float64[])
    for b5 in b5_values
        p = base_params(b5 = b5, a10 = a10, a11 = a11, b6 = b6, ctl_mode = ctl_mode)
        u0 = tumor_free_state(p; tumor_cells = 2.0)
        prob = SDEProblem(model_B!, sigma_B!, u0, (0.0, 365.0), p)

        success = 0
        tumor_means = Float64[]
        extinction_times = Float64[]

        for rep in 1:nreps
            sol = solve(prob, SOSRI(); saveat = 1.0, callback = stop_if_negative)
            tumor_traj = sol[1,:]
            extinct_idx = findfirst(x -> floor(x) < 1.0, tumor_traj)
            if isnothing(extinct_idx)
                success += 1
                push!(tumor_means, mean(tumor_traj))
            else
                push!(extinction_times, sol.t[extinct_idx])
            end
        end

        push!(out, (
            b5,
            success / nreps,
            isempty(tumor_means) ? NaN : mean(tumor_means),
            isempty(extinction_times) ? NaN : mean(extinction_times)
        ))
    end
    CSV.write(outfile, out)
    return out
end

function plot_b5_sde(df)
    p1 = plot(df.b5, df.success_prob, marker = :circle, lw = 3, xscale = :log10,
              xlabel = "b5", ylabel = "Success probability", title = "SDE metastasis establishment vs b5")
    p2 = plot(df.b5, df.mean_success_tumor, marker = :circle, lw = 3, xscale = :log10,
              xlabel = "b5", ylabel = "Mean tumor size (successful runs)", title = "Successful-run tumor burden vs b5")
    p3 = plot(df.b5, df.mean_time_to_extinction, marker = :circle, lw = 3, xscale = :log10,
              xlabel = "b5", ylabel = "Mean time to extinction", title = "Extinction time vs b5")
    return plot(p1, p2, p3, layout = (3,1), size = (900, 1000))
end


# main 
b5_vals = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3]
a10 = 2.0

traj = plot_b5_trajectories(b5_vals; a10 = a10, a11 = 0.0, b6 = 0.0,
                                  ctl_mode = :tumor_gated, init_mode = :tumor_free)
savefig(traj, "beta5_ode_trajectories.png")

ode_df = sweep_b5_ode(b5_vals; a10 = a10, a11 = 0.0, b6 = 0.0,
                      ctl_mode = :tumor_gated, init_mode = :tumor_free,
                      outfile = "beta5_ode_sweep.csv")
ode_plot = plot_b5_summary(ode_df)
savefig(ode_plot, "beta5_ode_summary.png")

sde_df = stochastic_b5_summary(b5_vals; a10 = a10, a11 = 0.0, b6 = 0.0,
                               ctl_mode = :tumor_gated, nreps = 100,
                               outfile = "beta5_sde_summary.csv")
sde_plot = plot_b5_sde(sde_df)
savefig(sde_plot, "beta5_sde_summary.png")
