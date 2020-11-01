#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["nacs_201804/data/data_20180409_194201.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "../../../experiments", iname))
               for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [(-80:8.0:120, 60:8.0:260)]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs)

const data_cs0 = datas_cs[1][1] # light off
const data_cs1 = datas_cs[1][2] # light on

const prefix = joinpath(@__DIR__, "../figures/pa_vectorshift_spectrum")

function rabiLine(det, t, Omega)
    Omega2 = Omega^2
    OmegaG2 = det^2 + Omega2
    return Omega2 / OmegaG2 * sin(√(OmegaG2) * t / 2)^2
end
function gen_rabi(t)
    return (x, p) -> p[1] .* rabiLine.(2π .* (x .- p[2]), t, p[3])
end
fit_cs0 = fit_survival(gen_rabi(17e-3), data_cs0, [0.9, 26, 100])
fit_cs1 = fit_survival(gen_rabi(17e-3), data_cs1, [0.9, 161, 100])
@show fit_cs0.uncs
@show fit_cs1.uncs

figure()
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_cs0.param[2], data_cs0),
                            fmt="C0o", label="PA off")
plot(fit_cs0.plotx .- fit_cs0.param[2], fit_cs0.ploty, "C0")
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v - fit_cs0.param[2], data_cs1),
                            fmt="C1o", label="PA on")
plot(fit_cs1.plotx .- fit_cs0.param[2], fit_cs1.ploty, "C1")
grid()
ylim([0, 1])
ylabel("F=3 Population")
xlabel("Detuning from Resonance (kHz)")
legend(loc="upper left", handlelength=0.6, handletextpad=0.3)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
