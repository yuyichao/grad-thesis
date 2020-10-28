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
        param, counts = matopen(joinpath(expdir, "pa_201909/data", "data_$name.mat")) do fd
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
    # datas = datas[[1:3:size(datas, 1) 2:3:size(datas, 1) 3:3:size(datas, 1)]]
    datas = datas[[1:2:size(datas, 1) 2:2:size(datas, 1)]]
    datas = NaCsData.map_params((i1, i2, v)->v[1], datas)
    return datas
end

const names_files = ["pa_201909/data/names_20190901_232500.mat"]
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

const prefix = joinpath(@__DIR__, "../figures/pa_spectrum_v12_blue_10")

fit_10 = fit_survival(model_lorentzian, datas[1][:, 1], [0.7, 0.35, -20, 50])
# fit_14 = fit_survival(model_lorentzian, datas[1][:, 2], [0.7, 0.35, -40, 50])

figure()
NaCsPlot.plot_survival_data(datas[1][:, 1], fmt="C0.", label="10 mW")
plot(fit_10.plotx, fit_10.ploty, "C0")
text(-155, 0.61, "\$\\Gamma_{10}=$(fit_10.uncs[4])\$ MHz", color="C0", fontsize="small")
grid()
ylim([0.28, 0.67])
# xlim([495.94, 496.16])
xlabel("PA Frequency (306496XXX MHz)")
ylabel("Two-Body Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
