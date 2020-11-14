#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using MAT
using PyPlot
using Statistics

const iname = joinpath(@__DIR__, "../../../experiments/nacs_atoms/data",
                       "data_20161019_171310.mat")

imgs, single_atom = matopen(iname) do fd
    read(fd, "Images"), read(fd, "SingleAtom") .!= 0
end

const whitelist = [27, 55, 78, 93, 144, 154, 171, 245, 311, 329, 542, 622, 629, 631, 632, 678,
                   897, 917, 918, 933, 961, 962, 1067, 1101, 1160, 1237, 1393, 1430, 1465,
                   1468, 1486, 1490, 1645, 1695, 1729, 1719]

const prefix = joinpath(@__DIR__, "../figures/loading_single_atoms")

figure()
implot = imshow(imgs[6:36, 6:36, 1695], interpolation="nearest",
                extent=[-5, 5, -5, 5])
implot[:set_cmap]("viridis")
xlabel("\$x\$ (\$\\mathrm{\\mu m}\$)")
ylabel("\$y\$ (\$\\mathrm{\\mu m}\$)")
annotate("Na", xy=(-0.5 / 1.5, 0), xytext=(-3.5, 3.5),
         arrowprops=Dict("color"=>"white", "shrink"=>0.1,
                         "width"=>6, "headwidth"=>12, "headlength"=>20),
         color="white", size=40, weight="bold")
annotate("Cs", xy=(2 / 1.5, 0.5 / 1.5), xytext=(2.5, -4),
         arrowprops=Dict("color"=>"white", "shrink"=>0.1,
                         "width"=>6, "headwidth"=>12, "headlength"=>20),
         color="white", size=40, weight="bold")
# colorbar()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
