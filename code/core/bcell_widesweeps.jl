using DifferentialEquations
using Plots
using Colors
using Statistics
using DataFrames
using CSV
using Random


# baseline parameters
function base_params(; a10 = 2.0, b5 = 1.0e-4, a11 = 1.0e-1, b6 = 1.0e-7,
                     b1 = 3.5e-6, b2 = 1.1e-7, ctl_mode = :direct)
    return (
        # original Kregerparameters
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

# initial states
function tumor_free_state(p; tumor_cells = 2.0)
    xM_star = p.a2/p.z2
    xN_star = p.z2*p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC_star = 0.0
    xB_star = p.a8/(p.z5 + p.b5 * xM_star)
    return [tumor_cells, xM_star, xN_star, xC_star, xB_star]
end

# model equations
function model_B!(du, u, p, t)
    xT, xM, xN, xC, xB = u

    B_to_CTL = if p.ctl_mode == :direct
        p.a11 * (xB / (p.g5 + xB))
    else
        p.a11 * (xT^2 / (p.g3 + xT^2)) * (xB / (p.g5 + xB))
    end

    # Tumor
    du[1] = p.a1 * xT * log(max(p.h / xT, 1.0)) - p.b1 * xT * xN - p.b2 * xT * xC - p.b6 * xT * xB - p.z1 * xT

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
        du[i] = u[i]
    end
end

function negative_condition(u, t, integrator)
    any(x -> x < 0.0, u)
end
stop_if_negative = DiscreteCallback(negative_condition, terminate!)

function solve_ode_model(; a10 = 2.0, b5 = 1e-4, a11 = 1e-1, b6 = 1e-7,
                         b1 = 3.5e-6, b2 = 1.1e-7, ctl_mode =:direct,
                         tmax = 365.0, saveat = 1.0)
    p = base_params(a10 = a10, b5 = b5, a11 = a11, b6 = b6, b1 = b1, b2 = b2, ctl_mode = ctl_mode)
    u0 = tumor_free_state(p; tumor_cells = 2.0)
    prob = ODEProblem(model_B!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5(); saveat = saveat, reltol = 1e-6, abstol = 1e-9,
                callback = stop_if_negative,
                isoutofdomain = (u,p,t)->any(x->x<0,u))
    return sol, p
end


# helper metrics
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

function run_one_parameter_sweep(param_name, values; base = (a10=2.0, b5=1e-4, a11=1e-1, b6=1e-7, b1=3.5e-6, b2=1.1e-7), ctl_mode=:direct)
    rows = DataFrame(param = Float64[], T60=Float64[], T100=Float64[], T365=Float64[], M365=Float64[], NK365=Float64[], CTL365=Float64[], B365=Float64[], t_to_1e5=Float64[], t_to_1e6=Float64[], auc_0_120=Float64[])
    for v in values
        kwargs = Dict(:a10 => base.a10, :b5 => base.b5, :a11 => base.a11, :b6 => base.b6, :b1 => base.b1, :b2 => base.b2, :ctl_mode => ctl_mode)
        kwargs[Symbol(param_name)] = v
        sol, _ = solve_ode_model(; kwargs...)
        m = metrics_from_sol(sol)
        push!(rows, (v, m.T60, m.T100, m.T365, m.M365, m.NK365, m.CTL365, m.B365, m.t_to_1e5, m.t_to_1e6, m.auc_0_120))
    end
    rename!(rows, :param => Symbol(param_name))
    return rows
end

function plot_one_parameter_summary(df, pname; baseline_value=nothing)
    x = df[!, Symbol(pname)]

    p1 = plot(x, df.T100, marker=:circle, lw=2.5, xscale=:log10,
              xlabel=pname, ylabel="Tumor at day 100", title="Tumor at day 100 vs $(pname)")
    p2 = plot(x, df.T365, marker=:circle, lw=2.5, xscale=:log10,
              xlabel=pname, ylabel="Tumor at day 365", title="Tumor at day 365 vs $(pname)")
    p3 = plot(x, df.t_to_1e5, marker=:circle, lw=2.5, xscale=:log10,
              xlabel=pname, ylabel="Time to 1e5 tumor cells", title="Time to threshold vs $(pname)")
    p4 = plot(x, df.auc_0_120, marker=:circle, lw=2.5, xscale=:log10,
              xlabel=pname, ylabel="Tumor AUC (0-120 d)", title="Early tumor burden vs $(pname)")
    return plot(p1, p2, p3, p4, layout=(2,2), size=(1000,800))
end

function heatmap_two_parameters(p1_name, p1_vals, p2_name, p2_vals;
                                base = (a10=2.0, b5=1e-4, a11=1e-1, b6=1e-7, b1=3.5e-6, b2=1.1e-7),
                                metric = :T100, ctl_mode=:direct)
    Z = zeros(length(p2_vals), length(p1_vals))
    for (i, p2v) in enumerate(p2_vals)
        for (j, p1v) in enumerate(p1_vals)
            kwargs = Dict(:a10 => base.a10, :b5 => base.b5, :a11 => base.a11, :b6 => base.b6, :b1 => base.b1, :b2 => base.b2, :ctl_mode => ctl_mode)
            kwargs[Symbol(p1_name)] = p1v
            kwargs[Symbol(p2_name)] = p2v
            sol, _ = solve_ode_model(; kwargs...)
            m = metrics_from_sol(sol)
            Z[i,j] = getproperty(m, metric)
        end
    end
    return heatmap(p1_vals, p2_vals, Z, xscale=:log10, yscale=:log10,
                   xlabel=p1_name, ylabel=p2_name,
                   title="$(String(metric)) for $(p1_name) × $(p2_name)", colorbar_title=String(metric)), Z
end

function run_small_sde_grid(pairs; nreps = 200, ctl_mode=:direct)
    rows = DataFrame(a10=Float64[], b5=Float64[], a11=Float64[], b6=Float64[], b1=Float64[], b2=Float64[], success_prob=Float64[])
    for par in pairs
        p = base_params(a10=par.a10, b5=par.b5, a11=par.a11, b6=par.b6, b1=par.b1, b2=par.b2, ctl_mode=ctl_mode)
        u0 = tumor_free_state(p; tumor_cells = 2.0)
        prob = SDEProblem(model_B!, sigma_B!, u0, (0.0, 365.0), p)
        success = 0
        for rep in 1:nreps
            sol = solve(prob, SOSRI(); saveat = 1.0, callback = stop_if_negative)
            tumor_traj = sol[1,:]
            extinct_idx = findfirst(x -> floor(x) < 1.0, tumor_traj)
            if isnothing(extinct_idx)
                success += 1
            end
        end
        push!(rows, (par.a10, par.b5, par.a11, par.b6, par.b1, par.b2, success / nreps))
    end
    return rows
end

# wider ranges
a10_vals = [1e-2, 3e-2, 1e-1, 3e-1, 1.0, 3.0, 10.0, 30.0, 100.0]
b5_vals  = [1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
a11_vals = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0, 3.0, 10.0]
b6_vals  = [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5]
b1_vals  = [1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4]
b2_vals  = [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5]

baseline = (a10=2.0, b5=1e-4, a11=1e-1, b6=1e-7, b1=3.5e-6, b2=1.1e-7)

# one-parameter wider sweeps
println(" wider one-parameter ODE sweeps")
for (pname, vals) in [("a10", a10_vals), ("b5", b5_vals), ("a11", a11_vals), ("b6", b6_vals), ("b1", b1_vals), ("b2", b2_vals)]
    df = run_one_parameter_sweep(pname, vals; base=baseline, ctl_mode=:direct)
    CSV.write("wide_sweep_$(pname).csv", df)
    fig = plot_one_parameter_summary(df, pname)
    savefig(fig, "wide_sweep_$(pname).png")
end

println("Running two-parameter heatmaps")
h1, _ = heatmap_two_parameters("a10", a10_vals, "b5", b5_vals; base=baseline, metric=:T100, ctl_mode=:direct)
savefig(h1, "heatmap_T100_a10_b5.png")

h2, _ = heatmap_two_parameters("a11", a11_vals, "b5", b5_vals; base=baseline, metric=:T100, ctl_mode=:direct)
savefig(h2, "heatmap_T100_a11_b5.png")

h3, _ = heatmap_two_parameters("a10", a10_vals, "b1", b1_vals; base=baseline, metric=:T100, ctl_mode=:direct)
savefig(h3, "heatmap_T100_a10_b1.png")

h4, _ = heatmap_two_parameters("a11", a11_vals, "b2", b2_vals; base=baseline, metric=:T100, ctl_mode=:direct)
savefig(h4, "heatmap_T100_a11_b2.png")