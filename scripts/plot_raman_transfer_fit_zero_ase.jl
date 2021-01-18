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
    if length(p) > 2
        p[1] .* exp.(.- x .* p[2]) .+ p[3]
    else
        p[1] .* exp.(.- x .* p[2])
    end
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
            return model_exp(t, gen_exp_param(p, etyp, idx, idx0))
        end
    end
    return wrapper.(xs)
end

module OneASE

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

import ..expdir, ..select_datas, ..model_exp, ..model_2d

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

fit = fit_survival(model, data_fit, [609.0, 2π * 1.5, 0, 2π / 0.2, 0, 0, 0.2,
                                     0.30, 0.8, 0.3, 0.3, 0.3],
                   plotx=false, lower=[-Inf; zeros(11)])
const param_1 = get_model_param(fit.param, 1)
const uncs_1 = get_model_param(fit.uncs, 1)
const scale_1 = 1 / (param_1[1] + param_1[2])

end

fit = fit_survival(model, data_fit, [594.29, 2π * 1.5, 0, 2π / 0.2, 0, 0, 0.1, 0.02,
                                     0.3, 0.3, 0.8, 0.3, 0.8, 0.3, 0.3, 0.3, 0.1], plotx=false)
const param_1 = get_model_param(fit.param, 1)
const uncs_1 = get_model_param(fit.uncs, 1)
const param_2 = get_model_param(fit.param, 2)
const scale_1 = 1 / (param_1[1] + param_1[2])
const scale_2 = 1 / (param_2[1] + param_2[2])

@show uncs_1 # 3.63(34)\times10^{-2}, 0.2824(52), 594.17(13), 25.28(37), 0.1490(93), 11.77(68)

const plot_time = linspace(0, 0.45, 1000)

figure()
NaCsPlot.plot_survival_data([OneASE.data_nacs_00; OneASE.data_nacs_t],
                            OneASE.scale_1, fmt="C0o")
plot(plot_time, model_2d.(plot_time, 609.0, (OneASE.param_1,)) .* OneASE.scale_1, "C0")
NaCsPlot.plot_survival_data([data_nacs_00; data_nacs_t], scale_1, fmt="C1o")
plot(plot_time, model_2d.(plot_time, 594.29, (param_1,)) .* scale_1, "C1")
errorbar([], [], [], fmt="C0o-", label="One Filter")
errorbar([], [], [], fmt="C1o-", label="Zero Filter")
xlim([0, 0.46])
legend(fontsize="small")
grid()
xlabel("Raman time (ms)")
ylabel("Fraction of Atomic Ground State")
NaCsPlot.maybe_save("$(prefix)_time")

NaCsPlot.maybe_show()
