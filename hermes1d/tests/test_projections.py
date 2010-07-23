from numpy import empty, pi, arange
from numpy.linalg import solve

from hermes1d.h1d_wrapper.h1d_wrapper import (assemble_projection_matrix_rhs,
        Mesh, FESolution)
from hermes_common._hermes_common import CooMatrix

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
    assemble_projection_matrix_rhs(m, A, rhs, projection_type="L2")
    x = solve(A.to_scipy_coo().todense(), rhs)
    sol_l2 = FESolution(m, x).to_discrete_function()
    A = CooMatrix(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, projection_type="H1")
    x = solve(A.to_scipy_coo().todense(), rhs)
    sol_h1 = FESolution(m, x).to_discrete_function()
    sol_l2.plot(False)
    sol_h1.plot(False)
