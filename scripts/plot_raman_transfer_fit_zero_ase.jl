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

const inames = ["nacs_202003/data/data_20200326_101325.mat",
                "nacs_202003/data/data_20200326_223438.mat",
                "nacs_202003/data/data_20200327_004910.mat",
                "nacs_202003/data/data_20200330_033004.mat",
                "nacs_202003/data/data_20200327_133537.mat",
                "nacs_202003/data/data_20200330_201614.mat",
                "nacs_202003/data/data_20200331_020659.mat",
                "nacs_202004/data/data_20200401_060025.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(expdir, iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [(593.0 .+ [-20; -6:2:6; 20], # 15 mW, 0.09 ms
                522.7 .+ [-15; -4.5:1.5:4.5; 15], # 12 mW, 0.12 ms
                447.5 .+ [-10; -3:0.6:3; 5], # 9 mW, 0.16 ms
                369.5 .+ [-5; -2:0.5:2; 5], # 6 mW, 0.27 ms
                287.5 .+ [-5; -0.6:0.15:0.6; 5], # 3 mW, 0.9 ms
                ),
               (593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.1 ms
                593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.14 ms
                ),
               ([0], # 15 mW, 0 ms
                593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.1 ms
                593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.14 ms
                ),
               (593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.32 ms
                369.5 .+ [-10; -2.5:0.5:2.5; 10], # 6 mW, 0.96 ms
                ),
               [0; [0.08, 0.16, 0.24, 0.32, 0.4, 0.48] .- 0.04], # 15 mW, 770.59429 MHz
               [0, 0.08, 0.16, 0.24, 0.32, 0.4, 0.48] .* 50,
               ([0.0, 2, 4, 8, 12, 20],
                [0.0, 2, 4, 8, 12, 20],
                [0.0, 20, 40, 80, 160],
                [0.0, 20, 40, 80, 160],
                [0.0, 20, 40, 80, 160] .* 2,
                [0.0, 20, 40, 80, 160] .* 2),
               ([0],
                [0],
                [0, 0.08, 0.16, 0.24, 0.3],
                )
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

const data_nacs_00 = datas_nacs[3][1]
const data_nacs_05 = datas_nacs[1][1]
const data_nacs_06 = [datas_nacs[2][1]; datas_nacs[3][2]]
const data_nacs_10 = [datas_nacs[2][2]; datas_nacs[3][3]]
const data_nacs_28 = datas_nacs[4][1]
const data_nacs_t = datas_nacs[5]

const data_pa_na1 = datas_na[6]
const data_pa_m1 = datas_nacs[6]

const data_pa_na2 = [datas_na[7][1]; datas_na[7][2]]
const data_pa_cs2 = datas_cs[7][1]
const data_pa_a2 = datas_nacs[7][1]
const data_pa_m2 = datas_nacs[7][2]

const data_nacs_m = datas_nacs[8][3]

const data_fit = [NaCsData.map_params((i, v) -> (1, v, 0.0, 1), data_nacs_00);
                  NaCsData.map_params((i, v) -> (1, v, 0.05, 2), data_nacs_05);
                  NaCsData.map_params((i, v) -> (1, v, 0.06, 1), data_nacs_06);
                  NaCsData.map_params((i, v) -> (1, v, 0.10, 1), data_nacs_10);
                  NaCsData.map_params((i, v) -> (1, v, 0.28, 1), data_nacs_28);
                  NaCsData.map_params((i, v) -> (1, 594.29, v, 1), data_nacs_t);
                  NaCsData.map_params((i, v) -> (2, v, 1, 3), data_pa_na1);
                  NaCsData.map_params((i, v) -> (2, v, 4, 4), data_pa_m1);
                  NaCsData.map_params((i, v) -> (2, v, 1, 5), data_pa_na2);
                  NaCsData.map_params((i, v) -> (2, v, 2, 6), data_pa_cs2);
                  NaCsData.map_params((i, v) -> (2, v, 3, 7), data_pa_a2);
                  NaCsData.map_params((i, v) -> (2, v, 4, 8), data_pa_m2);
                  NaCsData.map_params((i, v) -> (3, v, 5, 9, 2), data_nacs_m);
                  ]

const prefix = joinpath(@__DIR__, "../figures/raman_transfer_fit_zero_ase")

function model_exp(x, p)
    p[1] .* exp.(.- x .* p[2])
end

function model_expoff(x, p)
    p[3] .+ p[1] .* exp.(.- x .* p[2])
end

function get_model_param(p, idx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r = p
    p1 = p[8 + idx]
    p0 = p1 * p0r
    return (p0, p1, f0, Ω, Γ1 + Γ_na + Γ_cs, Γ2)
end

function gen_exp_param(p, typ, idx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs = p
    p1 = p[8 + idx]
    if typ == 1
        return (p1, Γ_na)
    elseif typ == 2
        return (p1, Γ_cs)
    elseif typ == 3
        return (p1, Γ_na + Γ_cs)
    elseif typ == 4
        return (p1, Γ1 + Γ_na + Γ_cs)
    elseif typ == 5
        return (p1, Γ2)
    end
end

function gen_exp_param(p, typ, idx, idx0)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs = p
    p0 = p[6 + idx0]
    p1 = p[8 + idx]
    if typ == 1
        return (p1, Γ_na, p0)
    elseif typ == 2
        return (p1, Γ_cs, p0)
    elseif typ == 3
        return (p1, Γ_na + Γ_cs, p0)
    elseif typ == 4
        return (p1, Γ1 + Γ_na + Γ_cs, p0)
    elseif typ == 5
        return (p1, Γ2, p0)
    end
end

function model(xs, p)
    # p: p0r, f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p1s...
    function wrapper(x)
        typ = x[1]
        if typ == 1
            _, f, t, idx = x
            return model_2d(t, f, get_model_param(p, idx))
        elseif typ == 2
            t = x[2]
            _, t, etyp, idx = x
            return model_exp(t, gen_exp_param(p, etyp, idx))
        else
            _, t, etyp, idx, idx0 = x
            return model_expoff(t, gen_exp_param(p, etyp, idx, idx0))
        end
    end
    return wrapper.(xs)
end
fit = fit_survival(model, data_fit, [594.29, 2π * 1.5, 0, 2π / 0.2, 0, 0, 0.1, 0.02,
                                     0.3, 0.3, 0.8, 0.3, 0.8, 0.3, 0.3, 0.3, 0.1], plotx=false)
const param_1 = get_model_param(fit.param, 1)
const uncs_1 = get_model_param(fit.uncs, 1)
const param_2 = get_model_param(fit.param, 2)
const scale_1 = 1 / (param_1[1] + param_1[2])
const scale_2 = 1 / (param_2[1] + param_2[2])

@show uncs_1 # 3.63(34)\times10^{-2}, 0.2824(52), 594.17(13), 25.28(37), 0.1490(93), 11.77(68)

ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, param_1)

const plot_freq_lo = 594.29 - 30
const plot_freq_hi = 594.29 + 30
const plot_freq = linspace(plot_freq_lo, plot_freq_hi, 1000)
const plot_time = linspace(0, 0.45, 1000)
const plot_pa_time = linspace(0, 25, 1000)

figure()
NaCsPlot.plot_survival_data([data_nacs_00; data_nacs_t], scale_1, fmt="C0o")
plot(plot_time, model_2d.(plot_time, 594.29, (param_1,)) .* scale_1, "C0")
# legend(fontsize="x-small", loc="upper right")
# text(0.09, 0.12, ("\$f_{res}=$(uncs_1[3] / 1000 + 770)\$ MHz\n" *
#                   "\$\\Omega_{Raman}=2\\pi\\times$(uncs_1[4] / 2π)\$ kHz\n" *
#                   "\$\\Gamma_{atom}=2\\pi\\times$(uncs_1[5] / 2π * 1000)\$ Hz\n" *
#                   "\$\\Gamma_{molecule}=2\\pi\\times$(uncs_1[6] / 2π)\$ kHz\n"),
#      color="C0", fontsize="small")
# annotate("Atom", (0, 0.32), xytext=(0.1, 0.34), arrowprops=Dict(:color=>"C3"), color="C3")
# annotate("Molecule", (0.147, 0.0344),
#          xytext=(0.08, 0.15), arrowprops=Dict(:color=>"C3"), color="C3")
# annotate("Atom", (0.26, 0.1),
#          xytext=(0.29, 0.17), arrowprops=Dict(:color=>"C3"), color="C3")
xlim([0, 0.46])
grid()
xlabel("Raman time (ms)")
ylabel("Fraction of Atomic Ground State")
NaCsPlot.maybe_save("$(prefix)_time")

NaCsPlot.maybe_show()