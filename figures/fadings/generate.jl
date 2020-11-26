#!/usr/bin/julia

using PyCall
using PyPlot

function draw_circle_fade(alpha_cb, X, Y, color, color2=color)
    img = Array{Float64}(undef, 4, X, Y)
    Xc = (1 + X) / 2
    Yc = (1 + Y) / 2
    Xr = Xc - 1
    Yr = Yc - 1
    for x in 1:X
        rx = (x - Xc) / Xr
        for y in 1:Y
            ry = (y - Yc) / Yr
            r = sqrt(rx^2 + ry^2)
            alpha = r <= 1 ? alpha_cb(r) : 0.0
            if alpha == 0
                img[:, x, y] .= 0
            else
                img[1, x, y], img[2, x, y], img[3, x, y] = color .* (1 - r) .+ color2 .* r
                img[4, x, y] = alpha
            end
        end
    end
    return img
end
function draw_linear_fade(alpha_cb, X, Y, color, color2=color)
    img = Array{Float64}(undef, 4, X, Y)
    for x in 1:X
        r = x / (X - 1)
        alpha = r <= 1 ? alpha_cb(r) : 0.0
        if alpha == 0
            img[:, x, :] .= 0
        else
            for y in 1:Y
                img[1, x, y], img[2, x, y], img[3, x, y] = color .* (1 - r) .+ color2 .* r
                img[4, x, y] = alpha
            end
        end
    end
    return img
end

save_img(name, array) = imsave(joinpath(@__DIR__, "$name.png"), PyReverseDims(array))

struct LinGrad
    r0::Float64
end

function (g::LinGrad)(r)
    if r <= g.r0
        return 1.0
    elseif r >= 1
        return 0.0
    else
        return (1 - r) / (1 - g.r0)
    end
end

struct PwrGrad
    pwr::Float64
end
function (g::PwrGrad)(r)
    if r <= 0
        return 1.0
    elseif r >= 1
        return 0.0
    else
        return 1 - r^g.pwr
    end
end

const grad_full = LinGrad(0)
const grad_glow = LinGrad(0.55)
const grad_glow2 = LinGrad(0.25)
const grad_fade_img = PwrGrad(3)

save_img("aom", draw_circle_fade(grad_full, 1000, 500, (1.0, 0.5, 0.0)))
save_img("lens", draw_circle_fade(grad_full, 600, 84, (0.0, 0.0, 1.0), (0.0, 0.6, 1.0)))
save_img("waveplate", draw_circle_fade(grad_full, 600, 42, (0.0, 0.0, 0.502), (0.0, 0.0, 0.8)))
save_img("pbs", draw_linear_fade(grad_full, 1000, 1, (0.0, 0.4, 1.0), (0.0, 0.502, 1.0)))
save_img("fade_img", draw_linear_fade(grad_fade_img, 1000, 1, (1.0, 1.0, 1.0)))
