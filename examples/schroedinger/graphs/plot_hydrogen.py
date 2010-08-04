#! /usr/bin/env python

from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel

import hydrogen_uniformpfem
import hydrogen_romanowski
import hydrogen_hpfem
import hydrogen_pfem
import hydrogen_pfem_uniform_init

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
            n, l, "k", "uniform $p$-FEM (L)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_romanowski.R_x[l], hydrogen_romanowski.R_y[n, l],
            n, l, "y", "$h$-FEM (Romanowski, L)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_pfem.R_x[l], hydrogen_pfem.R_y[n, l],
            n, l, "m", "$p$-FEM (L)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_pfem_uniform_init.R_x[l],
            hydrogen_pfem_uniform_init.R_y[n, l],
            n, l, "b", "$p$-FEM (U)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_hpfem.R_x[l], hydrogen_hpfem.R_y[n, l],
            n, l, "r", "$hp$-FEM (U)")
savefig("conv_l_0.png")
