#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
# using DataStructures
# using LsqFit
using LibArchive

using DelimitedFiles

const fname = joinpath(@__DIR__, "../data/NaCsPotentials.csv.zst")
const data = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    readdlm(reader, ',', skipstart=1)
end
const prefix = joinpath(@__DIR__, "../figures/pa_pes")

const wntofreq = 0.0299792458

figure()
axhline(0, linestyle="--", color="grey")
plot(data[:, 1], data[:, 2] .* wntofreq, "C0", label="\$X^1\\Sigma\$")
plot(data[:, 1], data[:, 3] .* wntofreq, "C1", label="\$a^3\\Sigma\$")
plot(data[:, 1], data[:, 4] .* wntofreq, "C2", label="\$A^1\\Sigma\$")
plot(data[:, 1], data[:, 5] .* wntofreq, "C3", label="\$b^3\\Pi\$")
plot(data[:, 1], data[:, 8] .* wntofreq, "C4", label="\$B^1\\Pi\$")
plot(data[:, 1], data[:, 9] .* wntofreq, "C5", label="\$c^3\\Sigma\$")
xlim([2.5, 9.9])
ylim([-200, 590])
text(5, -141, "\$\\mathrm{X^1\\Sigma}\$", color="C0", fontsize=15)
text(4.5, 16, "\$\\mathrm{a^3\\Sigma}\$", color="C1", fontsize=15)
text(6.95, 207, "\$\\mathrm{A^1\\Sigma}\$", color="C2", fontsize=15)
annotate("\$\\mathrm{b^3\\Pi}\$", (5.129, 234.7), xytext=(5.796, 229.3),
         arrowprops=Dict(:color=>"C3", :width=>0.1, :headwidth=>2, :headlength=>7),
         color="C3", fontsize=15)
text(4.4, 329.2, "\$\\mathrm{B^1\\Pi}\$", color="C4", fontsize=15)
text(3.9, 247.0, "\$\\mathrm{c^3\\Sigma}\$", color="C5", fontsize=15)
ax = gca()
ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)
xticks([3, 4, 5, 6, 7, 8, 9])
xlabel("Internuclear distance (\$\\mathrm{\\AA}\$)")
ylabel("Energy (THz)")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
