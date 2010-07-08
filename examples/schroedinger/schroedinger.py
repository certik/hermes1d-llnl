#! /usr/bin/env python

import sys
sys.path.insert(0, "../..")

from numpy import arange

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
        els.append((n, error))
    els.sort(key=lambda x: x[1])
    els.reverse()
    n, error = els[0]
    return n

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



def main():
    pts = arange(0, R, float(R)/(N_elem+1))
    orders = [P_init]*N_elem
    mesh = Mesh(pts, orders)
    for i in range(1):
        print "-"*80
        print "adaptivity iteration:", i
        N_dof = mesh.assign_dofs()
        A = CooMatrix(N_dof)
        B = CooMatrix(N_dof)
        assemble_schroedinger(mesh, A, B, l=0)
        eigs = solve_eig_numpy(A.to_scipy_coo(), B.to_scipy_coo())[:4]
        els2refine = []
        for E, eig in eigs:
            s = FESolution(mesh, eig)
            id = find_element_romanowski(s.get_element_coeffs())
            els2refine.append(id)
            print E, id
            #f = s.to_discrete_function()
            #print "plotting"
            #f.plot(False)
        els2refine = list(set(els2refine))
        print "Will refine the elements:", els2refine
        #plot_eigs(mesh, eigs)
        mesh = refine_mesh(mesh, els2refine)

if __name__ == "__main__":
    main()
