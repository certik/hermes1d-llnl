#! /usr/bin/env python

import sys
sys.path.insert(0, "../..")

from numpy import arange
from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel

from hermes1d import Mesh
from hermes1d.solvers.eigen import solve_eig_numpy, solve_eig_pysparse
from hermes1d.h1d_wrapper.h1d_wrapper import FESolution
from hermes1d.fekete.fekete import Function, Mesh1D
from hermes_common._hermes_common import CooMatrix

from _forms import assemble_schroedinger
from plot import plot_eigs, plot_file


N_elem = 4                         # number of elements
R = 150                            # right hand side of the domain
P_init = 6                         # initial polynomal degree

def find_element_romanowski(coeffs):
    els = []
    for n, e in enumerate(coeffs):
        error = min(abs(e[2:]))
        #print n, e, error
        els.append((n, error))
    els.sort(key=lambda x: x[1])
    els.reverse()
    n, error = els[0]
    return n, error

def refine_mesh(mesh, els2refine):
    new_pts = []
    pts, orders = mesh.get_mesh_data()
    new_pts.append(pts[0])
    for id in range(len(orders)):
        if id in els2refine:
            new_pts.append((pts[id]+pts[id+1])/2.)
        new_pts.append(pts[id+1])
    # assumes uniform order:
    orders = [orders[0]] * (len(new_pts)-1)
    return Mesh(new_pts, orders)

def do_plot(x, y, n, l):
    n_r = n - l - 1
    styles = {0: "-s", 1: "--o", 2: ":^", 3: "-.v"}
    plot(x, y, "k" + styles[n_r], label="$R_{%d%d}$" % (n, l))

    grid(True)
    ax = gca()
    xlabel("DOFs")
    ylabel("$E_{num}-E$")
    ax.set_yscale("log")
    title("l=%d" % l)
    legend()

def plot_conv(conv_graph, exact=None, l=None):
    assert exact is not None
    assert l is not None
    x = []
    y = [[], [], [], []]
    for dofs, energies in conv_graph:
        x.append(dofs)
        for i in range(4):
            y[i].append(energies[i]-exact[i])
    print x
    for i in range(4):
        n = l+1+i
        print n, l
        print y[i]
        do_plot(x, y[i], n, l)
    savefig("conv_l_0.png")

def main():
    #do_plot([23, 29, 41, 47], [0.1, 0.01, 0.001, 0.004], 1, 0)
    pts = arange(0, R, float(R)/(N_elem))
    pts = list(pts) + [R]
    #pts = list(pts) + [10000]
    orders = [P_init]*(len(pts)-1)
    mesh = Mesh(pts, orders)
    conv_graph = []
    l = 0
    error_tol = 1e-6
    #error_tol = 1e-2
    for i in range(100000):
        print "-"*80
        print "adaptivity iteration:", i
        N_dof = mesh.assign_dofs()
        pts, orders = mesh.get_mesh_data()
        print "Current mesh:"
        print pts
        A = CooMatrix(N_dof)
        B = CooMatrix(N_dof)
        assemble_schroedinger(mesh, A, B, l=l)
        eigs = solve_eig_numpy(A.to_scipy_coo(), B.to_scipy_coo())[:4]
        print
        els2refine = []
        errors = []
        energies = []
        for E, eig in eigs:
            s = FESolution(mesh, eig)
            id, error = find_element_romanowski(s.get_element_coeffs())
            els2refine.append(id)
            errors.append(error)
            energies.append(E)
            print E, id, error
            #f = s.to_discrete_function()
            #print "plotting"
            #f.plot(False)
        conv_graph.append((N_dof, energies))
        total_error = max(errors)
        print "Total error:", total_error
        if total_error < error_tol:
            break
        els2refine = list(set(els2refine))
        print "Will refine the elements:", els2refine
        mesh = refine_mesh(mesh, els2refine)
    plot_conv(conv_graph, exact=[-1./(2*n**2) for n in range(1+l, 5+l)], l=l)
    #plot_eigs(mesh, eigs)

if __name__ == "__main__":
    main()
