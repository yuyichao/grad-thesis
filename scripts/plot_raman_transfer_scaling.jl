#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival, get_plot_range
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const res_freqs = [668.35 21 772.1367 0.0082
                   668.35 21 772.1376 0.0049
                   668.35 21 772.1172 0.0048
                   668.35 18 771.924 0.033
                   668.35 12 771.323 0.050

                   560 15 770.59417 0.00013
                   560 6 770.369081 0.000037
                   560 3 770.287276 0.000015

                   565 15 770.61656 0.00021
                   565 6 770.378274 0.000080
                   565 3 770.291891 0.000038

                   555 15 770.57744 0.00013
                   555 6 770.362165 0.000043
                   555 3 770.283665 0.000017

                   503 15 770.483553 0.000049
                   503 12 770.43248 0.00011
                   503 9 770.379505 0.000028
                   503 6 770.324463 0.000019
                   503 3 770.2648307 0.0000084

                   492 15 770.47571 0.00010
                   492 12 770.42599 0.00032
                   492 9 770.375007 0.000069
                   492 6 770.321355 0.000037
                   492 3 770.263402 0.000028

                   510 15 770.50117 0.00083
                   510 12 770.44659 0.00021
                   510 9 770.390205 0.000079
                   510 6 770.331510 0.000057
                   510 3 770.258511 0.000023

                   605 15 770.76208 0.00027
                   605 12 770.65624 0.00046
                   605 9 770.54816 0.00058
                   605 6 770.437052 0.000072
                   605 3 770.321367 0.000038

                   625 15 770.89410 0.00024
                   625 6 770.490130 0.000059
                   625 3 770.347812 0.000022
                   ]

const rabi_freqs = [560 15 3.971 0.089
                    560 6 1.182 0.033
                    560 3 0.488 0.015

                    503 15 2.073 0.051
                    503 9 1.074 0.036
                    503 6 0.680 0.022
                    503 3 0.2544 0.0077

                    605 15 6.39 0.15
                    605 6 1.985 0.043
                    605 3 0.845 0.042

                    625 15 8.35 0.25
                    625 6 2.513 0.036
                    625 3 0.996 0.029
                    ]

function model_sqr(x, p)
    return p[1] .+ (p[2] .+ p[3] .* x) .* x
end

function model_omega(x, p)
    return p[1] .* x.^1.29
end

broacast_array(f, x) = f(x)
broacast_array(f, x::AbstractArray) = f.(x)

# x: freq, power, idx
# p: framan0, fpa0, offset, strength, bs...
function model_res(xs, p)
    function real_model(x)
        freq, power, idx = x
        framan0, fpa0, offset, strength = p
        b = idx == 0 ? 0.0 : p[4 + idx]
        p_sqr = (framan0, offset - strength / (freq - fpa0), b)
        return model_sqr(power, p_sqr)
    end
    return broacast_array(real_model, xs)
end

# x: freq, power
# p: offset, strength
function gen_model_rabi(fpa0)
    function model(xs, p)
        function real_model(x)
            freq, power = x
            offset, strength = p
            p_omega = (offset - strength / (freq - fpa0),)
            return model_omega(power, p_omega)
        end
        return broacast_array(real_model, xs)
    end
end

function gen_data(data_in)
    freq_ids = Dict{Float64,Int}()
    xs = Tuple{Float64,Float64,Int}[]
    ys = Float64[]
    uncs = Float64[]
    idmax = 0
    for i in 1:size(data_in, 1)
        freq = data_in[i, 1]
        power = data_in[i, 2]
        freq_id = get!(freq_ids, freq) do
            idmax += 1
            return idmax
        end
        freq_id = 1
        push!(xs, (freq, power, freq_id))
        push!(ys, data_in[i, 3])
        push!(uncs, data_in[i, 4])
    end
    return (x=xs, y=ys, unc=uncs, freqs=collect(keys(freq_ids)))
end

function fit_freq_model_res(model, data, _p0)
    p0 = [_p0; 0]
    return fit_data(model, data.x, data.y, data.unc, p0, plotx=false)
end

function get_plot_data_freq(data, fit, model, freq)
    local freq_id
    xs = Float64[]
    ys = Float64[]
    uncs = Float64[]
    for i in 1:length(data.x)
        x = data.x[i]
        x[1] == freq || continue
        if !@isdefined(freq_id)
            freq_id = x[3]
        else
            @assert(freq_id == x[3])
        end
        push!(xs, x[2])
        push!(ys, data.y[i])
        push!(uncs, data.unc[i])
    end
    function plot_func(x)
        return model((freq, x, freq_id), fit.param)
    end
    plotx = get_plot_range(xs)
    return (x=xs, y=ys, unc=uncs, plotx=plotx, ploty=plot_func.(plotx))
end

function get_plot_data_power(data, fit, model, power)
    xs = Float64[]
    xids = Float64[]
    ys = Float64[]
    uncs = Float64[]
    for i in 1:length(data.x)
        x = data.x[i]
        x[2] == power || continue
        push!(xs, x[1])
        push!(ys, data.y[i])
        push!(uncs, data.unc[i])
    end
    function plot_func(x)
        return model((x, power, 1), fit.param)
    end
    plotx = get_plot_range(xs)
    return (x=xs, y=ys, unc=uncs, plotx=plotx, ploty=plot_func.(plotx))
end

const data_res = gen_data(res_freqs)
const fit_res = fit_freq_model_res(model_res, data_res, [770.202, 705, 0, 100])
# 770.200516(24), 711.696(73), 2.194(17)\times10^{-3}, 4.1083(46), -1.9913(58)\times10^{-4}
@show fit_res.uncs

const model_rabi = gen_model_rabi(fit_res.param[2])
const data_rabi = gen_data(rabi_freqs)
const fit_rabi = fit_data(model_rabi, data_rabi.x, data_rabi.y, [0.1, 0.1], plotx=false)
# -6.83(27)\times10^{-2}, 28.00(33)
@show fit_rabi.uncs
const prefix = joinpath(@__DIR__, "../figures/raman_transfer_scaling")

const powers_plot = [15, 6, 3]

figure()
for i in 1:length(powers_plot)
    power = powers_plot[i]
    pd = get_plot_data_power(data_res, fit_res, model_res, power)
    errorbar(1000 ./ (pd.x .- fit_res.param[2]), (pd.y .- fit_res.param[1]) .* 1000,
             pd.unc .* 1000, fmt="C$(i - 1)o", label="$(power) mW")
    plot(1000 ./ (pd.plotx .- fit_res.param[2]), (pd.ploty .- fit_res.param[1]) .* 1000,
         "C$(i - 1)")
end
# p: framan0, fpa0, offset, strength, bs...
# text(560, 1700, ("\$f_{Raman0} + a\\cdot P\$\n" *
#                  "  \$-\\dfrac{b}{f-f_{PA0}}\\cdot P\$"))
# text(480, 740, ("\$f_{Raman0}=$(fit_res.uncs[1])\$ MHz\n" *
#                 "\$f_{PA0}=$(fit_res.uncs[2] + 288000)\$ GHz\n" *
#                 "\$a=$(fit_res.uncs[3] * 1000)\$ kHz/mW\n" *
#                 "\$b=$(fit_res.uncs[4])\$ MHz\$\\cdot\$GHz/mW"), fontsize="small")
legend(loc="upper right", fontsize="small", handlelength=0.6, handletextpad=0.3)
grid()
yticks([0, 200, 400, 600])
ylim([0, 790])
xlabel("\$1 / \\Delta~(\\mathrm{THz^{-1}})\$")
ylabel("Light Shift (kHz)")
NaCsPlot.maybe_save("$(prefix)_light_shift_det")

figure()
for i in 1:length(powers_plot)
    power = powers_plot[i]
    pd = get_plot_data_power(data_rabi, fit_rabi, model_rabi, power)
    plot(1000 ./ (pd.plotx .- fit_res.param[2]), pd.ploty, "C$(i - 1)")
    errorbar(1000 ./ (pd.x .- fit_res.param[2]), pd.y, pd.unc,
             fmt="C$(i - 1)o", label="$(power) mW")
end
# p: offset, strength
# text(500, 5.0, ("\$\\left(a-\\dfrac{b}{f-705 GHz}\\right)\\cdot P^{1.29}\$"))
# text(500, 2.8, ("\$a=$(fit_rabi.uncs[1] * 1000)\$ Hz/mW\$^{1.29}\$\n" *
#               "\$b=$(fit_rabi.uncs[2])\$ kHz\$\\cdot\$GHz/mW\$^{1.29}\$"), fontsize="small")
legend(loc="upper right", fontsize="small", handlelength=0.6, handletextpad=0.3)
grid()
yticks([0, 2, 4, 6, 8])
# xticks([-200, -160, -120, -80])
ylim([0, 9.7])
xlabel("\$1 / \\Delta~(\\mathrm{THz^{-1}})\$")
ylabel("\$\\Omega_{\\mathrm{R}}~(2\\pi\\cdot \\mathrm{kHz})\$")
NaCsPlot.maybe_save("$(prefix)_rabi_det")

const freq = 560

figure()
plot2 = get_plot_data_freq(data_res, fit_res, model_res, freq)
errorbar(plot2.x[1], (plot2.y[1] .- fit_res.param[1]) .* 1000,
         plot2.unc[1] .* 1000, fmt="C0o")
errorbar(plot2.x[2], (plot2.y[2] .- fit_res.param[1]) .* 1000,
         plot2.unc[2] .* 1000, fmt="C1o")
errorbar(plot2.x[3], (plot2.y[3] .- fit_res.param[1]) .* 1000,
         plot2.unc[3] .* 1000, fmt="C2o")
plot(plot2.plotx, (plot2.ploty .- fit_res.param[1]) .* 1000, "k")
grid()
xlabel("Tweezer Power (mW)")
ylabel("Light Shift (kHz)")
xlim([2.4, 17])
ylim([77, 440])
NaCsPlot.maybe_save("$(prefix)_light_shift_560")

figure()
ax2 = gca()
plot2 = get_plot_data_freq(data_rabi, fit_rabi, model_rabi, freq)
errorbar(plot2.x[1], plot2.y[1], plot2.unc[1], fmt="C0o")
errorbar(plot2.x[2], plot2.y[2], plot2.unc[2], fmt="C1o")
errorbar(plot2.x[3], plot2.y[3], plot2.unc[3], fmt="C2o")
plot(plot2.plotx, plot2.ploty, "k")
xscale("log")
yscale("log")
grid()
xticks([3, 6, 12], ["3", "6", "12"])
yticks([0.5, 1, 2, 4], ["0.5", "1", "2", "4"])
ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax2.set_xticks(3:17, minor=true)
xlabel("Tweezer Power (mW)")
ylabel("\$\\Omega_{\\mathrm{R}} (2\\pi\\cdot \\mathrm{kHz})\$")
xlim([2.4, 17])
ylim([0.42, 4.2])
NaCsPlot.maybe_save("$(prefix)_rabi_560")

NaCsPlot.maybe_show()
