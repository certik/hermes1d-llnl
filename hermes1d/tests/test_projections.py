from numpy import empty, pi, arange, array, sin, cos
from numpy.linalg import solve

from hermes1d.h1d_wrapper.h1d_wrapper import (assemble_projection_matrix_rhs,
        Mesh, FESolution, Function as Function2)
from hermes1d.fekete.fekete import Function, Mesh1D
from hermes_common._hermes_common import CooMatrix

class FunctionSin(Function2):

    def eval_f(self, x):
        return sin(x)

    def eval_dfdx(self, x):
        return cos(x)

f_sin = FunctionSin()


def test_l2_h1_proj_run():
    """
    Test that the projections run.

    It doesn't test if it's correct.
    """
    pts = arange(0, 2*pi, 1)
    orders = [3]*(len(pts)-1)
    m = Mesh(pts, orders)
    n_dof = m.assign_dofs()
    A = CooMatrix(n_dof)
    rhs = empty(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="L2")
    x = solve(A.to_scipy_coo().todense(), rhs)
    sol_l2 = FESolution(m, x).to_discrete_function()
    A = CooMatrix(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="H1")
    x = solve(A.to_scipy_coo().todense(), rhs)
    sol_h1 = FESolution(m, x).to_discrete_function()
    sol_l2.plot(False)
    sol_h1.plot(False)

def test_l2_h1_proj1():
    """
    Tests the correctness of the projections.
    """
    pts = arange(0, 2*pi, 3)
    orders = [2]*(len(pts)-1)
    m = Mesh(pts, orders)

    pts = array(list(arange(0, pts[-1], 0.1)) + [pts[-1]])
    orders = [6]*(len(pts)-1)
    f_exact = Function(lambda x: sin(x), Mesh1D(pts, orders))

    n_dof = m.assign_dofs()
    A = CooMatrix(n_dof)
    rhs = empty(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="L2")
    x = solve(A.to_scipy_coo().todense(), rhs)
    sol_l2 = FESolution(m, x).to_discrete_function()
    A = CooMatrix(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="H1")
    x = solve(A.to_scipy_coo().todense(), rhs)
    sol_h1 = FESolution(m, x).to_discrete_function()
    assert (sol_l2 - f_exact).l2_norm() < 0.07
    assert (sol_h1 - f_exact).l2_norm() < 0.07

def test_l2_h1_proj2():
    """
    Tests the correctness of the projections.
    """
    pts = arange(0, 2*pi, 3)
    orders = [4]*(len(pts)-1)
    m = Mesh(pts, orders)

    pts = array(list(arange(0, pts[-1], 0.1)) + [pts[-1]])
    orders = [6]*(len(pts)-1)
    f_exact = Function(lambda x: sin(x), Mesh1D(pts, orders))

    n_dof = m.assign_dofs()
    A = CooMatrix(n_dof)
    rhs = empty(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="L2")
    x = solve(A.to_scipy_coo().todense(), rhs)
    sol_l2 = FESolution(m, x).to_discrete_function()
    A = CooMatrix(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="H1")
    x = solve(A.to_scipy_coo().todense(), rhs)
    sol_h1 = FESolution(m, x).to_discrete_function()
    assert (sol_l2 - f_exact).l2_norm() < 0.002
    assert (sol_h1 - f_exact).l2_norm() < 0.002

def test_l2_h1_proj3():
    """
    Tests conversion to FE basis.
    """
    pts = arange(0, 2*pi, 0.1)
    orders = [2]*(len(pts)-1)
    m = Mesh(pts, orders)

    f = Function(lambda x: sin(x), Mesh1D(pts, orders))

    n_dof = m.assign_dofs()
    A = CooMatrix(n_dof)
    rhs = empty(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f, projection_type="L2")
    x = solve(A.to_scipy_coo().todense(), rhs)
    sol_l2 = FESolution(m, x).to_discrete_function()
    A = CooMatrix(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f, projection_type="H1")
    x = solve(A.to_scipy_coo().todense(), rhs)
    sol_h1 = FESolution(m, x).to_discrete_function()
    assert sol_l2 == f
    assert sol_h1 == f
