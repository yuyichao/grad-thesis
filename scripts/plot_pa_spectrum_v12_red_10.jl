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
            pl = read(fd, "ParamList")
            sg = read(fd, "ScanGroup")
            f = sg["base"]["params"]["fWavemeter0"]
            amps = sg["base"]["vars"]["params"]["Merge"]["Cs"]["Wait"]["PTwzr"]
            c = read(fd, "SingleAtomLogical") .!= 0
            return [(f, amps[Int(pl[i])]) for i in 1:size(c, 3)], c
        end
        data = NaCsData.select_count(param, counts, selector)
        if !@isdefined datas
            datas = data
        else
            datas = [datas; data]
        end
    end
    datas = datas[[1:2:size(datas, 1) 2:2:size(datas, 1)]]
    datas = NaCsData.map_params((i1, i2, v)->v[1], datas)
    return datas
end

const names_files = ["pa_201908/data/names_20190830_005224.mat"]
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
    ds[1] = NaCsData.map_params((i1, i2, v)->(v - 496) * 1000, ds[1])
end

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end

const prefix = joinpath(@__DIR__, "../figures/pa_spectrum_v12_red_10")

# fit_5 = fit_survival(model_lorentzian, datas[1][:, 1], [0.7, 0.35, 60, 60])
fit_10 = fit_survival(model_lorentzian, datas[1][:, 2], [0.7, 0.35, 100, 167], plot_lo=-20)

figure()
NaCsPlot.plot_survival_data(datas[1][:, 2], fmt="C0.")
plot(fit_10.plotx, fit_10.ploty, "C0-")
text(10, 0.71, "\$\\Gamma=$(fit_10.uncs[4])\$ MHz", color="C0")
grid()
ylim([0.4, 0.8])
yticks([0.4, 0.5, 0.6, 0.7, 0.8])
# xlim([495.9, 496.34])
xlabel("PA Frequency (306496XXX MHz)")
ylabel("Two-Body Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
