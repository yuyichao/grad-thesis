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
        param, counts = matopen(joinpath(expdir, "nacs_201911/data", "data_$name.mat")) do fd
            pl = read(fd, "ParamList")
            sg = read(fd, "ScanGroup")
            f = sg["base"]["params"]["fWavemeter1"]
            # ts = sg["scans"]["vars"]["params"]["TMergeWait"]
            # f = pl[1]
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

const names_files = ["nacs_201911/data/names_20191124_122128.mat"]
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
    # ds[1] = [ds[1]; ds[2]]
    # resize!(ds, 1)
    ds[1] = NaCsData.map_params((i, v)->v - 307000, ds[1])
end

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end

const prefix = joinpath(@__DIR__, "../figures/pa_spectrum_v13")

# fit = fit_survival(model_lorentzian, datas[1], [0.8, .7, 492.35, 0.5])

figure()
const ptrprops = Dict(:width=>2, :headlength=>6, :headwidth=>6, :color=>"C1")
# axvline(879.75, color="C1", ls="--")
NaCsPlot.plot_survival_data(datas[1], fmt="C0.-")
annotate("", xy=(879.75, 0.16), xytext=(879.75, 0.11), arrowprops=ptrprops)
grid()
ylim([0.10, 0.80])
xlabel("PA Frequency (307XXX GHz)")
ylabel("Two-Body Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
