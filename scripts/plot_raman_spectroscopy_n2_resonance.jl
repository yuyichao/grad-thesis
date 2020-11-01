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

const inames = ["nacs_201810/data/data_20181103_165031.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "../../../experiments", iname))
               for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [(81:0.06:90) .* 2 .- 60]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2,), (3, 4,)), maxcnts, specs)
const data_nacs = datas_nacs[1]

const prefix = joinpath(@__DIR__, "../figures/raman_spectroscopy_n2_resonance")

figure()
NaCsPlot.plot_survival_data(data_nacs, fmt="C0o-") # 10 mW
grid()
xlabel("Raman Frequency (MHz)")
ylabel("Atom Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
