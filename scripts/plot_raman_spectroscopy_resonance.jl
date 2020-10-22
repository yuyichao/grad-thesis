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

const inames = ["nacs_201810/data/data_20181114_230217.mat",
                "nacs_201810/data/data_20181116_113355.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "../../../experiments", iname))
               for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [(1.0:(10 * 41), 0.0:0.0),
               (1.0:(10 * 41), 0.0:0.0)]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2,), (3, 4,)), maxcnts, specs)
const data_nacs = [datas_nacs[1][1]; datas_nacs[2][1]]

const times = Float64[1, 2, 4, 6, 8, 12, 16, 20, 25, 30]
const freqs = 298.62:0.002:298.70

const freq_datas =
    [NaCsData.map_params((i, v)->freqs[i], data_nacs[i:10:end]) for i in 1:10]

const prefix = joinpath(@__DIR__, "../figures/raman_spectroscopy_resonance")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
fit2 = fit_survival(model_lorentzian, freq_datas[2], [0.8, 0.5, 298.665, 0.01])

figure()
# 2 ms
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->(v - fit2.param[3]) .* 1000,
                                                freq_datas[2]), fmt="C0o")
plot((fit2.plotx .- fit2.param[3]) .* 1000, fit2.ploty, "C0")
grid()
ylim([0.2, 0.8])
xlabel("Detuning from Resonance (kHz)")
ylabel("Atom Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
