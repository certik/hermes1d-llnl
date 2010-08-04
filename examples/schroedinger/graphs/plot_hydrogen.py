#! /usr/bin/env python

from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel

import hydrogen_uniformpfem

def do_plot(x, y, n, l, color="k", label=""):
    n_r = n - l - 1
    if n_r == 0:
        plot(x, y, color + "-o", label=label)
    else:
        plot(x, y, color + "-o")

    grid(True)
    ax = gca()
    xlabel("DOFs")
    ylabel("$E_{num}-E$")
    ax.set_yscale("log")
    title("l=%d" % l)
    legend()

n_eig = 3
l = 0
print "Saving to conv_l_0.png"
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_uniformpfem.R_x[l], hydrogen_uniformpfem.R_y[n, l],
            n, l, "k", "uniform $p$-FEM")
savefig("conv_l_0.png")
