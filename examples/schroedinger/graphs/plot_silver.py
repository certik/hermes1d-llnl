#! /usr/bin/env python

from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel

import silver_uniformpfem50
import silver_uniformpfem4
import silver_pfem
import silver_hpfem

def do_plot(x, y, n, l, color="k", label=""):
    n_r = n - l - 1
    if n_r == 0:
        plot(x, y, color + "-", label=label)
    else:
        plot(x, y, color + "-")

    grid(True)
    ax = gca()
    xlabel("DOFs")
    ylabel("$E_{num}-E$")
    ax.set_yscale("log")
    title("l=%d" % l)
    legend()

n_eig = 50
l = 0
print "Saving to conv_l_0.png"
for i in range(n_eig):
    n = l+1+i
    do_plot(silver_uniformpfem50.R_x[l], silver_uniformpfem50.R_y[n, l],
            n, l, "y", "uniform $p$-FEM, 50 elms")
for i in range(n_eig):
    n = l+1+i
    do_plot(silver_uniformpfem4.R_x[l], silver_uniformpfem4.R_y[n, l],
            n, l, "k", "uniform $p$-FEM, 4 elms")
for i in range(n_eig):
    n = l+1+i
    do_plot(silver_pfem.R_x[l], silver_pfem.R_y[n, l],
            n, l, "b", "$p$-FEM")
for i in range(n_eig):
    n = l+1+i
    do_plot(silver_hpfem.R_x[l], silver_hpfem.R_y[n, l],
            n, l, "r", "$hp$-FEM")
savefig("conv_l_0.png")
