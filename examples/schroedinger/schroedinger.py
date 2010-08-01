#! /usr/bin/env python

import sys
sys.path.insert(0, "../..")

from numpy import arange, empty, zeros, array
from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel

from hermes1d import Mesh
from hermes1d.solvers.eigen import solve_eig_numpy, solve_eig_pysparse, \
        solve_eig_scipy
from hermes1d.h1d_wrapper.h1d_wrapper import FESolution, calc_error_estimate, \
        calc_solution_norm, adapt
from hermes1d.fekete.fekete import Function, Mesh1D
from hermes_common._hermes_common import CooMatrix

from _forms import assemble_schroedinger
from plot import plot_eigs, plot_file


N_elem = 50                         # number of elements
R = 150                            # right hand side of the domain
P_init = 3                         # initial polynomal degree
error_tol = 1e-8                   # error tolerance
eqn_type="R"                      # either R or rR
NORM = 1 # 1 ... H1; 0 ... L2;
THRESHOLD = 0.7
#error_tol = 1e-2

def find_element_romanowski(coeffs):
    """
    Finds the smallest coefficient at each element (error) and return the
    element with the largest error.

    Effectivelly it is just using the coefficient at the highest bubble
    function and it needs at least quadratic elements (or higher) to work.
    """
    els = []
    for n, e in enumerate(coeffs):
        error = min(abs(e[2:]))
        #print n, e, error
        els.append((n, error))
    els.sort(key=lambda x: x[1])
    els.reverse()
    n, error = els[0]
    return n, error

def refine_mesh_romanowski(mesh, solutions):
    """
    Uses Romanowski refinement for all solutions in 'solutions'.

    Solutions are given as vectors coming from the matrix solver.
    """
    els2refine = []
    errors = []
    for sol in solutions:
        s = FESolution(mesh, sol)
        id, error = find_element_romanowski(s.get_element_coeffs())
        els2refine.append(id)
        errors.append(error)
    els2refine = list(set(els2refine))
    print "Will refine the elements:", els2refine
    mesh = refine_mesh(mesh, els2refine)
    return mesh

def refine_mesh_h1_adapt(mesh, solutions):
    """
    Uses H1 adaptivity refinement for all solutions in 'solutions'.

    Solutions are given as vectors coming from the matrix solver.
    """
    # so far only for one solution:
    assert len(solutions) == 1
    sol = solutions[0]
    mesh_ref = mesh.reference_refinement()
    print "Fine mesh created (%d DOF)." % mesh_ref.get_n_dof()
    return mesh_ref, [1.0]

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
    styles = {0: "-s", 1: "--o", 2: ":^", 3: "-.v", 4: "-.v", 5: "-.v"}
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
    n_eig = len(conv_graph[0][1])
    x = []
    y = [[] for n in range(n_eig)]
    for dofs, energies in conv_graph:
        x.append(dofs)
        for i in range(n_eig):
            y[i].append(energies[i]-exact[i])
    f = open("data.py", "w")
    f.write("R_x = {\n")
    f.write("        %d: %s,\n" % (l, x))
    f.write("    }\n")
    f.write("R_y = {\n")
    for i in range(n_eig):
        n = l+1+i
        f.write("        (%d, %d): %s,\n" % (n, l, y[i]))
        #do_plot(x, y[i], n, l)
    f.write("    }\n")
    #savefig("conv_l_0.png")

def flip_vectors(mesh, eigs, mesh_ref, eigs_ref, test_it=False):
    x_c = 1e-3
    for i in range(len(eigs)):
        s = FESolution(mesh, eigs[i])
        s_ref = FESolution(mesh_ref, eigs_ref[i])
        if s.value(x_c) < 0:
            #print "  Multiplying %d-th coarse eigenvector by (-1)" % i
            eigs[i] = -eigs[i]
        if s_ref.value(x_c) < 0:
            #print "  Multiplying %d-th ref. eigenvector by (-1)" % i
            eigs_ref[i] = -eigs_ref[i]

        if test_it:
            # Test it:
            s = FESolution(mesh, eigs[i]).to_discrete_function()
            s_ref = FESolution(mesh_ref, eigs_ref[i]).to_discrete_function()
            same_norm = (s-s_ref).l2_norm()
            flipped_norm = (s+s_ref).l2_norm()
            print same_norm, flipped_norm
            if same_norm > flipped_norm:
                c = min(same_norm, flipped_norm) / max(same_norm, flipped_norm)
                print "Warning: the flip is wrong, c=", c
                # If "c" is almost one, then the vectors can't really be
                # aligned anyway:
                assert c > 0.9
                #s.plot(False)
                #s_ref.plot()

def solve_schroedinger(mesh, l=0, Z=1, eqn_type=eqn_type, eig_num=4):
    """
    Solves the Schroedinger equation on the given mesh.

    Returns the energies and eigenfunctions.
    """
    # TODO: return the eigenfunctions as FESolutions
    N_dof = mesh.assign_dofs()
    A = CooMatrix(N_dof)
    B = CooMatrix(N_dof)
    assemble_schroedinger(mesh, A, B, l=l, Z=Z, eqn_type=eqn_type)
    eigs = solve_eig_scipy(A.to_scipy_coo(), B.to_scipy_coo())
    eigs = eigs[:eig_num]
    assert len(eigs) == eig_num
    energies = [E for E, eig in eigs]
    eigs = [eig for E, eig in eigs]
    return N_dof, array(energies), eigs

def adapt_mesh(mesh, eigs, l=0, Z=1, adapt_type="hp"):
    """
    Adapts the mesh using the adaptivity type 'adapt_type'.

    Returns a new instance of the H1D mesh.

    adapt_type .... one of: h, hp, p, uniform-p, romanowski
    """
    if adapt_type == "romanowski":
        m = refine_mesh_romanowski(mesh, eigs)
        pts, orders = m.get_mesh_data()
        return Mesh(pts, orders)
    elif adapt_type == "uniform-p":
        pts, orders = mesh.get_mesh_data()
        orders = array(orders) + 1
        return Mesh(pts, orders)
    elif adapt_type in ["h", "p", "hp"]:
        mesh_ref = mesh.reference_refinement()
        print "Fine mesh created (%d DOF)." % mesh_ref.get_n_dof()
        N_dof, energies, eigs_ref = solve_schroedinger(mesh_ref, l=l, Z=Z,
                eqn_type=eqn_type, eig_num=len(eigs))
        flip_vectors(mesh, eigs, mesh_ref, eigs_ref)
        print "    Done."
        sols = []
        sols_ref = []
        print "Normalizing solutions..."
        for i in range(len(eigs)):
            e = (eigs[i]).copy()
            coarse_h1_norm = FESolution(mesh, e).h1_norm()
            e /= coarse_h1_norm
            sols.append(e)
            e = (eigs_ref[i]).copy()
            reference_h1_norm = FESolution(mesh_ref, e).h1_norm()
            e /= reference_h1_norm
            sols_ref.append(e)
            #print "H1 norms:"
            #print "coarse    (%d):" % i, coarse_h1_norm
            #print "reference (%d):" % i, reference_h1_norm
        print "    Done."
        meshes = []
        mesh_orig = mesh.copy()
        mesh_orig.assign_dofs()
        errors = []
        for sol, sol_ref in zip(sols, sols_ref):
            mesh = mesh_orig.copy()
            mesh.assign_dofs()
            mesh_ref = mesh.reference_refinement()
            mesh_ref.assign_dofs()
            mesh.copy_vector_to_mesh(sol, 0)
            mesh_ref.copy_vector_to_mesh(sol_ref, 0)
            err_est_total, err_est_array = calc_error_estimate(NORM, mesh, mesh_ref)
            ref_sol_norm = calc_solution_norm(NORM, mesh_ref)
            err_est_rel = err_est_total/ref_sol_norm
            print "Relative error (est) = %g %%\n" % (100.*err_est_rel)
            errors.append(err_est_rel)
            # TODO: adapt using all the vectors:
            # 0 ... hp, 1 ... h, 2 ... p
            if adapt_type == "hp":
                ADAPT_TYPE = 0
            elif adapt_type == "h":
                ADAPT_TYPE = 1
            elif adapt_type == "p":
                ADAPT_TYPE = 2
            else:
                raise ValueError("Unkown adapt_type")
            adapt(NORM, ADAPT_TYPE, THRESHOLD, err_est_array, mesh, mesh_ref)
            meshes.append(mesh)
        pts, orders = mesh_orig.get_mesh_data()
        mesh = Mesh1D(pts, orders)
        for m in meshes:
            pts, orders = m.get_mesh_data()
            m = Mesh1D(pts, orders)
            mesh = mesh.union(m)
        pts, orders = mesh.get_mesh_data()
        mesh = Mesh(pts, orders)
        return mesh
    else:
        raise ValueError("Unknown adapt_type")


def main():
    #do_plot([23, 29, 41, 47], [0.1, 0.01, 0.001, 0.004], 1, 0)
    #pts = arange(0, R, float(R)/(N_elem))
    #pts = list(pts) + [R]
    par = 100.
    a, b, = 0., 150.
    Ne = N_elem
    r = par**(1./(Ne-1))
    pts = [(r**i-1)/(r**Ne-1)*(b-a)+a for i in range(Ne+1)]
    #pts = list(pts) + [10000]
    orders = [P_init]*(len(pts)-1)
    #pts = (0, 4.6875, 9.375, 18.75, 23.4375, 28.125, 32.8125, 35.15625, 37.5,
    #        42.1875, 46.875, 56.25, 65.625, 75.0, 150)
    #orders = (10, 8, 9, 4, 3, 3, 2, 2, 3, 2, 3, 2, 1, 1)
    #pts = (0, 4.6875, 9.375, 18.75, 23.4375, 28.125, 32.8125, 35.15625, 37.5,
    #        42.1875, 46.875, 49.21875, 51.5625, 56.25, 65.625, 75.0, 150)
    #orders = (10, 8, 9, 4, 4, 4, 2, 2, 3, 2, 2, 1, 2, 2, 1, 2)
    mesh = Mesh(pts, orders)
    conv_graph = []
    l=0
    Z = 47
    exact_energies=[-1.*Z**2/(2*n**2) for n in range(1+l,50+1+l)]
    old_energies = None
    for i in range(1000000):
        print "-"*80
        print "adaptivity iteration:", i
        if eqn_type == "rR":
            mesh.set_bc_left_dirichlet(0, 0)
            mesh.set_bc_right_dirichlet(0, 0)
        pts, orders = mesh.get_mesh_data()
        print "Current mesh:"
        print pts
        print orders
        N_dof, energies, eigs = solve_schroedinger(mesh, l=l, Z=Z,
                eqn_type=eqn_type, eig_num=50)
        conv_graph.append((N_dof, energies))
        # This doesn't work well:
        #if old_energies is not None:
        #    err = max(abs(old_energies - energies))
        #    if err < error_tol:
        #        print "Maximum error in energies:", err
        #        break
        err = max(energies - exact_energies)
        print "Maximum error in energies:", err
        if err < error_tol:
            break
        old_energies = energies
    #    exact_energies = array(exact_energies)
    #    print energies - exact_energies
        mesh = adapt_mesh(mesh, eigs, l=l, Z=Z, adapt_type="hp")
    plot_conv(conv_graph, exact=exact_energies, l=l)
    #plot_eigs(mesh, eigs)

if __name__ == "__main__":
    main()
