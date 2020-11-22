#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

include("../../../experiments/nacs_202003/molecular_raman_model.jl")

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const expdir = joinpath(@__DIR__, "../../../experiments")

const inames = ["nacs_202008/data/data_20200830_004406.mat",
                "nacs_202008/data/data_20200830_183835.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(expdir, iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [([0.0], # 15 mW, 0 ms
                609.0 .+ (-8:1.6:8), # 15 mW, 0.13 ms
                609.0 .+ (-12:2.4:12), # 15 mW, 0.26 ms
                [0.02, 0.04, 0.07, 0.085, 0.10, 0.13, 0.145, 0.16, 0.175, 0.19, 0.22,
                 0.25] .* 1.6 .- 0.01, # 15 mW, 770.609 MHz
                ),
               ([0.0], # 6 mW, 0 ms
                327.55 .+ (-2.5:0.5:2.5), # 6 mW, 0.51 ms
                327.55 .+ (-3:0.6:3), # 6 mW, 1.01 ms
                [0.04, 0.07, 0.085, 0.10, 0.115, 0.13, 0.16, 0.175, 0.19, 0.22, 0.25, 0.265,
                 0.28] .* 5.5 .- 0.01, # 6 mW, 770.32755 MHz
                [0.0, 2, 4, 10],
                [0.0, 20, 40, 80, 120],
                [0.0, 20, 40, 80, 160] .* 2,
                [0.0, 25, 50, 100],
                [0.0, 50, 100, 200],
                [0.0, 50, 100, 200] .* 2),
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((-1, 2), (-3, 4,)), maxcnts, specs)
const datas_na = select_datas(datas, NaCsData.select_single((1, -2), (3, -4,)), maxcnts, specs)

const data_nacs_t = datas_nacs[1][4] # Survival 1
const data_nacs_00 = datas_nacs[1][1] # Survival 1
const data_nacs_12 = datas_nacs[1][2] # Survival 1
const data_nacs_25 = datas_nacs[1][3] # Survival 1

const data_pa_na = [datas_na[2][4 + 4]; datas_na[2][1 + 4]]
const data_pa_cs = datas_cs[2][4 + 4]
const data_pa_a = datas_nacs[2][4 + 4]
const data_pa_m = datas_nacs[2][1 + 4]

const data_fit = [NaCsData.map_params((i, v) -> (1, 609.0, v, 1), data_nacs_00);
                  NaCsData.map_params((i, v) -> (1, v, 0.12, 1), data_nacs_12);
                  NaCsData.map_params((i, v) -> (1, v, 0.25, 1), data_nacs_25);
                  NaCsData.map_params((i, v) -> (1, 609.0, v, 1), data_nacs_t);
                  NaCsData.map_params((i, v) -> (2, v, 0.0, 2), data_pa_na);
                  NaCsData.map_params((i, v) -> (3, v, 0.0, 3), data_pa_cs);
                  NaCsData.map_params((i, v) -> (4, v, 0.0, 4), data_pa_a);
                  NaCsData.map_params((i, v) -> (5, v, 0.0, 5), data_pa_m);]

const prefix = joinpath(@__DIR__, "../figures/raman_transfer_fit_one_ase")

function model_exp(x, p)
    if length(p) > 2
        p[1] .* exp.(.- x .* p[2]) .+ p[3]
    else
        p[1] .* exp.(.- x .* p[2])
    end
end

function get_model_param(p, idx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r = p
    p1 = p[7 + idx]
    p0 = p1 * p0r
    return (p0, p1, f0, Ω, Γ1 + Γ_na + Γ_cs, Γ2)
end

function gen_pa_param(p, typ, pidx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r = p
    p1 = p[7 + pidx]
    if typ == 2
        return (p1, Γ_na)
    elseif typ == 3
        return (p1, Γ_cs)
    elseif typ == 4
        return (p1, Γ_na + Γ_cs)
    elseif typ == 5
        return (p1, Γ1 + Γ_na + Γ_cs)
    end
end

function model(xs, p)
    # p: f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r, p1s...
    function wrapper(x)
        typ = x[1]
        if typ == 1
            _, f, t, idx = x
            return model_2d(t, f, get_model_param(p, idx))
        else
            _, t, _, idx = x
            return model_exp(t, gen_pa_param(p, typ, idx))
        end
    end
    return wrapper.(xs)
end
function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end

fit = fit_survival(model, data_fit, [609.0, 2π * 1.5, 0, 2π / 0.2, 0, 0, 0.2,
                                     0.30, 0.8, 0.3, 0.3, 0.3],
                   plotx=false, lower=[-Inf; zeros(11)])
const param_1 = get_model_param(fit.param, 1)
const uncs_1 = get_model_param(fit.uncs, 1)
const scale_1 = 1 / (param_1[1] + param_1[2])
const param_pa = gen_pa_param(fit.param, 5, 5)

@show uncs_1 # 2.19(47)\times10^{-2}, 0.2684(80), 608.23(11), 24.95(56), 0.154(14), 8.19(65)

const data_nacs_12′ = NaCsData.map_params((i, v)->v - param_1[3], data_nacs_12)
const data_nacs_25′ = NaCsData.map_params((i, v)->v - param_1[3], data_nacs_25)

ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, param_1)

fit_l12 = let
    local xs, ratios, uncs = NaCsData.get_values(data_nacs_12)
    xs = copy(xs)
    ratios = ratios[:, 2]
    uncs = uncs[:, 2]
    push!(xs, param_1[3] - 200)
    push!(xs, param_1[3] + 200)
    push!(ratios, 1 / scale_1)
    push!(ratios, 1 / scale_1)
    push!(uncs, uncs_1[2].s)
    push!(uncs, uncs_1[2].s)
    fit_data(model_lorentzian, xs, ratios, uncs, [param_1[2], 0.20, 609, 5],
             plot_lo=param_1[3] - 15, plot_hi=param_1[3] + 15)
end
const data_nacs_12_shift = NaCsData.map_params((i, v)->v - fit_l12.param[3], data_nacs_12)

const plot_freq_lo = 609.0 - 13
const plot_freq_hi = 609.0 + 13
const plot_freq = linspace(plot_freq_lo, plot_freq_hi, 1000)
const plot_time = linspace(0, 0.4, 1000)
const plot_pa_time = linspace(0, 21, 1000)

figure()
NaCsPlot.plot_survival_data(data_nacs_12_shift, scale_1, fmt="C0.", label="0.12 ms")
plot(fit_l12.plotx .- fit_l12.param[3], fit_l12.ploty .* scale_1, "C0")
x₋ = fit_l12.param[3] - fit_l12.param[4] / 2
x₊ = fit_l12.param[3] + fit_l12.param[4] / 2
y_pm = (model_lorentzian(x₋, fit_l12.param) + model_lorentzian(x₊, fit_l12.param)) / 2 .* scale_1
ax = gca()
ax.annotate("", (-fit_l12.param[4] / 2 - 0.2, y_pm), xytext=(-fit_l12.param[4] / 2 - 5, y_pm),
            arrowprops=Dict(:color=>"C3"))
ax.annotate("", (fit_l12.param[4] / 2 + 0.2, y_pm), xytext=(fit_l12.param[4] / 2 + 5, y_pm),
            arrowprops=Dict(:color=>"C3"))
text(0, 0.12 .* scale_1, "\$\\mathbf{\\Gamma_{FWHM}=$(fit_l12.uncs[4]) kHz}\$",
     fontsize=17, color="C3", ha="center")
ylim([0.05, 1])
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("Detuning from resonance (kHz)")
ylabel("Fraction of Atomic Ground State")
NaCsPlot.maybe_save("$(prefix)_freq1_fwhm")

figure()
errorbar([0.0], [ratio_00 * scale_1], [uncs_00 * scale_1], fmt="C0o", label="0.00 ms")
plot([plot_freq_lo, plot_freq_hi] .- param_1[3], [v00 * scale_1, v00 * scale_1], "C0-")
NaCsPlot.plot_survival_data(data_nacs_12′, scale_1, fmt="C1o", label="0.12 ms")
plot(plot_freq .- param_1[3], model_2d.(0.12, plot_freq, (param_1,)) .* scale_1, "C1")
NaCsPlot.plot_survival_data(data_nacs_25′, scale_1, fmt="C2o", label="0.25 ms")
plot(plot_freq .- param_1[3], model_2d.(0.25, plot_freq, (param_1,)) .* scale_1, "C2")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("Raman Detuning (kHz)")
ylabel("Fraction of Atomic Ground State")
NaCsPlot.maybe_save("$(prefix)_freq")

figure()
NaCsPlot.plot_survival_data([data_nacs_00; data_nacs_t], scale_1, fmt="C0o")
plot(plot_time, model_2d.(plot_time, 609.0, (param_1,)) .* scale_1, "C0")
# text(0.10, 0.54, ("\$f_{res}=$(uncs_1[3] / 1000 + 770)\$ MHz\n" *
#                   "\$\\Omega_{Raman}=2\\pi\\times$(uncs_1[4] / 2π) kHz\$\n" *
#                   "\$\\Gamma_{atom}=2\\pi\\times$(uncs_1[5] / 2π * 1000)\$ Hz\n" *
#                   "\$\\Gamma_{molecule}=2\\pi\\times$(uncs_1[6] / 2π) kHz\$"),
#      color="C0", fontsize="small")
xlim([0, 0.4])
grid()
xlabel("Raman Time (ms)")
ylabel("Fraction of Atomic Ground State")
NaCsPlot.maybe_save("$(prefix)_time")

figure()
NaCsPlot.plot_survival_data(data_pa_m, 1 / param_pa[1], fmt="C0o")
plot(plot_pa_time, model_exp(plot_pa_time, param_pa) ./ param_pa[1], "C0-")
xlim([0, 17])
ylim([0, 1])
grid()
xlabel("Atom Hold Time (ms)")
ylabel("Two-Body Survival")
NaCsPlot.maybe_save("$(prefix)_atom_time")

NaCsPlot.maybe_show()
