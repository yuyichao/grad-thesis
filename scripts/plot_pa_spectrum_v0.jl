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
using MAT

const expdir = joinpath(@__DIR__, "../../../experiments")

function load_data(names, selector)
    local datas
    for name in names
        param, counts = matopen(joinpath(expdir, "nacs_201804/data", "data_$name.mat")) do fd
            pl = read(fd, "ParamList")
            # sg = read(fd, "ScanGroup")
            # f = sg["base"]["params"]["fWavemeter1"]
            # ts = sg["scans"]["vars"]["params"]["TMergeWait"]
            f = pl[1]
            c = read(fd, "SingleAtomLogical") .!= 0
            return [f for i in 1:size(c, 3)], c
        end
        data = NaCsData.select_count(param, counts, selector)
        if !@isdefined datas
            datas = data
        else
            datas = [datas; data]
        end
    end
    # datas = datas[[1:2:size(datas, 1) 2:2:size(datas, 1)]]
    # datas = NaCsData.map_params((i1, i2, v)->v[1], datas)
    return datas
end

const names_files = ["nacs_201804/data/names_20180417_234100.mat"
                     "nacs_201804/data/names_20180418_071502.mat"]
const datas = [load_data(matopen(fd->read(fd, "names"), joinpath(expdir, names_file)),
                         NaCsData.select_single((1, 2,), (3, 4,)))
               for names_file in names_files]
const datas_na = [load_data(matopen(fd->read(fd, "names"), joinpath(expdir, names_file)),
                            NaCsData.select_single((1, 2,), (3,)))
                  for names_file in names_files]
const datas_cs = [load_data(matopen(fd->read(fd, "names"), joinpath(expdir, names_file)),
                            NaCsData.select_single((1, 2,), (4,)))
                  for names_file in names_files]
const datas_nana = [load_data(matopen(fd->read(fd, "names"), joinpath(expdir, names_file)),
                              NaCsData.select_single((1, -2,), (3,)))
                    for names_file in names_files]
const datas_cscs = [load_data(matopen(fd->read(fd, "names"), joinpath(expdir, names_file)),
                              NaCsData.select_single((-1, 2,), (4,)))
                    for names_file in names_files]

for ds in [datas_cs, datas_na, datas, datas_cscs, datas_nana]
    ds[1] = NaCsData.map_params((i, v)->v * 11.576904 + 606.768489, ds[1])
    ds[2] = NaCsData.map_params((i, v)->v * 11.420860 + 608.195733, ds[2])
    # ds[1] = [ds[1]; ds[2]]
    # resize!(ds, 1)
end

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end

const prefix = joinpath(@__DIR__, "../figures/pa_spectrum_v0")

figure()
const ptrprops = Dict(:width=>2, :headlength=>6, :headwidth=>6, :color=>"C1")
# axvline(693.33, color="C1", ls="--")
# axvline(694.29, color="C1", ls="--")
# axvline(698.03, color="C1", ls="--")
# axvline(698.55, color="C1", ls="--")
# axvline(705.39, color="C1", ls="--")
NaCsPlot.plot_survival_data(datas[1], fmt="C0.-")
NaCsPlot.plot_survival_data(datas[2], fmt="C0.-")
annotate("", xy=(693.33, 0.17), xytext=(693.33, 0.12), arrowprops=ptrprops)
annotate("", xy=(694.29, 0.08), xytext=(694.29, 0.03), arrowprops=ptrprops)
annotate("", xy=(698.03, 0.15), xytext=(698.03, 0.10), arrowprops=ptrprops)
annotate("", xy=(698.55, 0.11), xytext=(698.55, 0.16), arrowprops=ptrprops)
annotate("", xy=(705.39, 0.08), xytext=(705.39, 0.03), arrowprops=ptrprops)
grid()
ylim([0, 0.95])
xlabel("PA Frequency (288XXX GHz)")
ylabel("Two-Body Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
