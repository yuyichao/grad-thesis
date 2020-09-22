#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsCalc
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

function fit(model, data, p0; plotx=nothing, use_unc=false, plot_scale=1.1)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
    end
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        plotx = linspace(mid - span * plot_scale / 2, mid + span * plot_scale / 2, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios[:, 2], 1 ./ uncs[:, 2].^2, p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
    return (param=fit.param, unc=estimate_errors(fit),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

model_peak(x, p) = p[1] .* exp.(.-((x .- p[2]) ./ p[3]).^2)

const data_dir = joinpath(@__DIR__, "..", "..", "..", "experiments", "na_rsc_201801", "data")

const iname_a = joinpath(data_dir, "data_20180121_145134.mat")
const iname_b = joinpath(data_dir, "data_20180122_214615.mat")
const iname_c = joinpath(data_dir, "data_20180123_155129.mat")

const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const data_a = NaCsData.select_count(params_a, logicals_a, selector)
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const data_b = NaCsData.select_count(params_b, logicals_b, selector)
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)
const data_c = NaCsData.select_count(params_c, logicals_c, selector)

# With cooling
const spec_a = OrderedDict(
    :x=>(linspace(18.11, 18.20, 10), linspace(19.05, 19.20, 16),
         linspace(19.54, 19.69, 16)),
    :y=>(linspace(18.10, 18.20, 11), linspace(19.04, 19.20, 17),
         linspace(19.52, 19.71, 20)),
    :z=>(linspace(18.49, 18.545, 12), linspace(18.57, 18.635, 14),
         linspace(18.66, 18.71, 11), linspace(18.75, 18.80, 11),
         linspace(18.83, 18.89, 13), linspace(18.91, 18.97, 13),
         linspace(18.99, 19.06, 15)),
    :x0=>linspace(18.50, 18.70, 11),
    :y0=>linspace(18.50, 18.70, 11),
)

# Without cooling
const spec_b = OrderedDict(
    :xp1=>linspace(0, 135, 16),
    :yp1=>linspace(6, 90, 15),
    :zp1=>linspace(17, 255, 15),
    :x0=>linspace(3.5, 52.5, 15),
    :y0=>linspace(2.5, 37.5, 15),
    :z0=>linspace(7, 105, 15),
    :xf=>linspace(18, 20, 101),
    :yf=>linspace(18, 20, 101),
    :zf=>linspace(18.5, 19.1, 121)
)

const spec_c = OrderedDict(
    :xp1=>linspace(0, 135, 16),
    :yp1=>linspace(6, 90, 15),
    :zp1=>linspace(17, 255, 15),
    :x0=>linspace(3.5, 52.5, 15),
    :y0=>linspace(2.5, 37.5, 15),
    :z0=>linspace(7, 105, 15),
    :xf=>linspace(18, 20, 101),
    :yf=>linspace(18, 20, 101),
    :zf=>linspace(18.4, 18.55, 31),
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)
const split_c = NaCsData.split_data(data_c, spec_c)

const prefix = joinpath(@__DIR__, "../figures/na_rsc_spectrum")

to_sideband(f) = (i, v)->(v - f) * 1000

const data_cold_rx = NaCsData.map_params(to_sideband(18.649), split_a[:x])
const data_cold_ry = NaCsData.map_params(to_sideband(18.648), split_a[:y])
const data_cold_az = NaCsData.map_params(to_sideband(18.605), split_a[:z])

const data_hot_rx = NaCsData.map_params(to_sideband(18.649), split_b[:xf])
const data_hot_ry = NaCsData.map_params(to_sideband(18.648), split_b[:yf])
const data_hot_az = NaCsData.map_params(to_sideband(18.598), [split_b[:zf][1:(end - 10)];
                                                              split_c[:zf][15:end]])

fig = figure()

ax0 = fig[:add_subplot](111)    # The big subplot
ax0[:spines]["top"][:set_color]("none")
ax0[:spines]["bottom"][:set_color]("none")
ax0[:spines]["left"][:set_color]("none")
ax0[:spines]["right"][:set_color]("none")
ax0[:tick_params](labelcolor="w", top="off", bottom="off", left="off", right="off")
ylabel("F=1 population")
ax0[:get_yaxis]()[:set_label_coords](-0.105, 0.5)

ax1 = fig[:add_subplot](211)
# Radial x
# Without cooling
NaCsPlot.plot_survival_data(data_hot_rx, fmt="C3o-", label="\$x\$ initial")
# With cooling
NaCsPlot.plot_survival_data(data_cold_rx[1], fmt="C0s-", label="\$x\$ cooled")
NaCsPlot.plot_survival_data(data_cold_rx[2], fmt="C0s-")
NaCsPlot.plot_survival_data(data_cold_rx[3], fmt="C0s-")
axhline(0.96, c="#99aaaa", ls="-.")
grid()
ylim([0, 1])
xlim([-700, 1400])
legend(fontsize=15, loc="center right")
setp(ax1[:get_xticklabels](), visible=false)
ax1[:tick_params](axis="x", length=0)
ax1[:get_yaxis]()[:set_label_coords](-0.105, 0.5)
yticks([0, 0.25, 0.5, 0.75, 1.0], ["", "0.25", "0.50", "0.75", "1.00"])

ax2 = fig[:add_subplot](212)
subplots_adjust(hspace=0.0)
# Without cooling
NaCsPlot.plot_survival_data(data_hot_ry, fmt="C3o-", label="\$y\$ initial")
# With cooling
NaCsPlot.plot_survival_data(data_cold_ry[1], fmt="C0s-", label="\$y\$ cooled")
NaCsPlot.plot_survival_data(data_cold_ry[2], fmt="C0s-")
NaCsPlot.plot_survival_data(data_cold_ry[3], fmt="C0s-")
axhline(0.96, c="#99aaaa", ls="-.")
grid()
ylim([0, 1])
xlim([-700, 1400])
legend(fontsize=15, loc="center right")
xlabel("\$\\delta\$, Detuning from carrier (kHz)")
ax2[:get_yaxis]()[:set_label_coords](-0.105, 0.5)
yticks([0, 0.25, 0.5, 0.75, 1.0], ["0.00", "0.25", "0.50", "0.75", "1.00"])

NaCsPlot.maybe_save("$(prefix)_r")

fig = figure()
# Without cooling
NaCsPlot.plot_survival_data(data_hot_az, fmt="C3o-", label="\$z\$ initial")

# With cooling
NaCsPlot.plot_survival_data(data_cold_az[1], fmt="C0s-", label="\$z\$ cooled")
NaCsPlot.plot_survival_data(data_cold_az[2], fmt="C0s-")
NaCsPlot.plot_survival_data(data_cold_az[3], fmt="C0s-")
NaCsPlot.plot_survival_data(data_cold_az[4], fmt="C0s-")
NaCsPlot.plot_survival_data(data_cold_az[5], fmt="C0s-")
NaCsPlot.plot_survival_data(data_cold_az[6], fmt="C0s-")
NaCsPlot.plot_survival_data(data_cold_az[7], fmt="C0s-")
axhline(0.96, c="#99aaaa", ls="-.")
grid()
ylim([0, 1])
xlim([-145, 470])
legend(fontsize=15, loc="center right")
xlabel("\$\\delta\$, Detuning from carrier (kHz)")
# xlabel("δ, Detuning from carrier (kHz)")
ylabel("F=1 population")
NaCsPlot.maybe_save("$(prefix)_az")

NaCsPlot.maybe_show()
