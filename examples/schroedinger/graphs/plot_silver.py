#! /usr/bin/env python

from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel

import silver_uniformpfem
import silver_pfem

def do_plot(x, y, n, l, color="k"):
    n_r = n - l - 1
    plot(x, y, color + "-", label="$R_{%d%d}$" % (n, l))

    grid(True)
    ax = gca()
    xlabel("DOFs")
    ylabel("$E_{num}-E$")
    ax.set_yscale("log")
    title("l=%d" % l)
    #legend()

n_eig = 50
l = 0
for i in range(n_eig):
    n = l+1+i
    do_plot(silver_uniformpfem.R_x[l], silver_uniformpfem.R_y[n, l], n, l)
    do_plot(silver_pfem.R_x[l], silver_pfem.R_y[n, l], n, l, "b")
savefig("conv_l_0.png")
