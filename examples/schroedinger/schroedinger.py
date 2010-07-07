#! /usr/bin/env python

from hermes1d import Mesh, CooMatrix

from _forms import assemble_schroedinger
from utils import solve_eig_numpy, solve_eig_pysparse
from plot import plot_eigs, plot_file

N_eq = 1
N_elem = 40                         # number of elements
R = 20                              # right hand side of the domain
P_init = 2                          # initial polynomal degree

def main():
    mesh = Mesh(0, R, N_elem, P_init, N_eq)
    N_dof = mesh.assign_dofs()
    A = CooMatrix(N_dof)
    B = CooMatrix(N_dof)
    assemble_schroedinger(mesh, A, B)
    eigs = solve_eig_numpy(A.to_scipy_coo(), B.to_scipy_coo())
    #eigs = solve_eig_pysparse(A.to_scipy_coo(), B.to_scipy_coo())
    E, v = eigs[0]
    plot_eigs(mesh, eigs)

if __name__ == "__main__":
    main()
