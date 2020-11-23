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

const fname = joinpath(@__DIR__, "../../damop-2020/20MHz_linewidth_c3SigmaOnly_3322_80kHzConfinement.csv.zst")
const data = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    readdlm(reader, ',', skipstart=1)
end
const prefix = joinpath(@__DIR__, "../figures/raman_transfer_theory_ratio")

figure()
ax1 = gca()
l1 = plot(data[:, 1] .- 288625.081, abs.(data[:, 2] ./ 2π / 1000), "C0",
          label="\$\\Omega_{\\mathrm{R}}\$")
l2 = plot(data[:, 1] .- 288625.081, abs.(data[:, 4] ./ 2π / 1000), "C1",
          label="\$\\Gamma_{\\mathrm{s}}\$")
text(-32, 2.5, "v'=0", fontsize="small")
ylabel("\$\\Omega_{\\mathrm{R}}\$ and \$\\Gamma_{\\mathrm{s}}\$ (\$2\\pi\\cdot \\mathrm{kHz})\$",
       fontsize="small")
xlabel("One-Photon Detuning (GHz)")
ylim([0, 27.5])
grid()
tax1 = ax1.twinx()
l3 = tax1.plot(data[:, 1] .- 288625.081, abs.(data[:, 2] ./ data[:, 4]), "C2",
               label="\$\\frac{\\Omega_{\\mathrm{R}}}{\\Gamma_{\\mathrm{s}}}\$", linewidth=2)
ls = [l1; l2; l3]
legend(ls, [l.get_label() for l in ls], fontsize="small",
       loc="lower right", bbox_to_anchor=(1, 0.05))
xlim([-35, 35])
ylim([0, 55])
ylabel("\$\\Omega_{\\mathrm{R}}/\\Gamma_{\\mathrm{s}}\$", fontsize="small", color="C2")
tax1.tick_params(axis="y", labelcolor="C2")
setp(tax1.get_yticklabels(), fontweight="bold")
NaCsPlot.maybe_save("$(prefix)_v0")

figure()
ax2 = gca()
l1 = plot(data[:, 1] .- 351271.53, abs.(data[:, 2] ./ 2π / 1000), "C0",
          label="\$\\Omega_{\\mathrm{R}}\$")
l2 = plot(data[:, 1] .- 351271.53, abs.(data[:, 4] ./ 2π / 1000), "C1",
          label="\$\\Gamma_{\\mathrm{s}}\$")
text(-32, 130, "v'=63", fontsize="small")
ylabel("\$\\Omega_{\\mathrm{R}}\$ and \$\\Gamma_{\\mathrm{s}}\$ (\$2\\pi\\cdot \\mathrm{kHz})\$",
       fontsize="small")
xlabel("One-Photon Detuning (GHz)")
ylim([0, 596])
grid()
tax2 = ax2.twinx()
l3 = tax2.plot(data[:, 1] .- 351271.53, abs.(data[:, 2] ./ data[:, 4]), "C2",
              label="\$\\frac{\\Omega_{\\mathrm{R}}}{\\Gamma_{\\mathrm{s}}}\$", linewidth=2)
ls = [l1; l2; l3]
legend(ls, [l.get_label() for l in ls], fontsize="small", loc="upper left")
xlim([-35, 35])
ylim([0, 1.49])
ylabel("\$\\Omega_{\\mathrm{R}}/\\Gamma_{\\mathrm{s}}\$", fontsize="small", color="C2")
tax2.tick_params(axis="y", labelcolor="C2")
setp(tax2.get_yticklabels(), fontweight="bold")
NaCsPlot.maybe_save("$(prefix)_vhi")

NaCsPlot.maybe_show()
