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

const pos = [0, -20, 20, 40, 60, 80, 100, 120] .* 0.140
const shift = [28.2, 9, 82, 110, 146, 136, 67, 30] # kHz
const shift_err = [1.3, 2, 2, 2, 2, 2, 2, 1.4] # kHz

const prefix = joinpath(@__DIR__, "../figures/pa_vectorshift_align")

function model_gaussian(x, p)
    p[1] ./ exp.(((x .- p[2]) ./ p[3]).^2)
end

fit = fit_data(model_gaussian, pos, shift, [140, 8.0, 6])
@show fit.uncs

figure()
errorbar(pos .- fit.param[2], shift, shift_err, fmt="C0o")
plot(fit.plotx .- fit.param[2], fit.ploty, "C0")
grid()
ylim([0, ylim()[2]])
ylabel("Resonance Shift")
xlabel("Distance from center (\$\\mathrm{\\mu m}\$)")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
