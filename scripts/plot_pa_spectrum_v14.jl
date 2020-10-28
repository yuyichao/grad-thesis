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
        param, counts = matopen(joinpath(expdir, "pa_201908/data", "data_$name.mat")) do fd
            sg = read(fd, "ScanGroup")
            f = sg["base"]["params"]["fWavemeter"]
            c = read(fd, "SingleAtomLogical") .!= 0
            return [f for _ in 1:size(c, 3)], c
        end
        data = NaCsData.select_count(param, counts, selector)
        if !@isdefined datas
            datas = data
        else
            datas = [datas; data]
        end
    end
    return datas
end

const names_files = ["pa_201908/data/names_20190816_204109.mat"]
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

const prefix = joinpath(@__DIR__, "../figures/pa_spectrum_v14")

for ds in [datas_cs, datas_na, datas, datas_cscs, datas_nana]
    ds[1] = NaCsData.map_params((i, v)->v - 3000, ds[1])
end

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
fit = fit_survival(model_lorentzian, datas[1], [0.8, 0.3, 255.45, 0.13])

figure()
const ptrprops = Dict(:width=>2, :headlength=>6, :headwidth=>6, :color=>"C1")
NaCsPlot.plot_survival_data(datas[1], fmt="C0.-")
annotate("", xy=(fit.param[3], 0.41), xytext=(fit.param[3], 0.38), arrowprops=ptrprops)
grid()
ylim([0.37, 0.75])
xlabel("PA Frequency (309XXX GHz)")
ylabel("Two-Body Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
