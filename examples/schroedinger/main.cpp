#include "hermes1d.h"

#include "python_api.h"
#include "_hermes1d_api.h"

#include "forms.h"

static int N_eq = 1;
int N_elem = 40;                         // number of elements
double A = 0, B = 20;                     // domain end points
int P_init = 2;                           // initial polynomal degree

/******************************************************************************/
int main(int argc, char* argv[]) {
  // create mesh
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  // you can set the zero dirichlet at the right hand side
  //mesh.set_bc_right_dirichlet(0, 0);

  // variable for the total number of DOF
  int N_dof = mesh->assign_dofs();

  CooMatrix *mat1 = new CooMatrix(N_dof);
  CooMatrix *mat2 = new CooMatrix(N_dof);
  int l=0;
  assemble_schroedinger(mesh, mat1, mat2, l);


  Python p;

  p.exec("print 'Python initialized'");
  p.push("A", c2py_CooMatrix(mat1));
  p.push("B", c2py_CooMatrix(mat2));
  p.exec("from utils import solve_eig_numpy, solve_eig_pysparse");
  p.exec("eigs = solve_eig_numpy(A.to_scipy_coo(), B.to_scipy_coo())");
  //p.exec("eigs = solve_eig_pysparse(A.to_scipy_coo(), B.to_scipy_coo())");
  p.exec("E, v = eigs[0]");

  if (import_hermes1d___hermes1d())
      throw std::runtime_error("Can't import hermes1d");
  p.push("mesh",  c2py_Mesh(mesh));
  printf("2\n");
  p.exec("from plot import plot_eigs, plot_file");
  p.exec("plot_eigs(mesh, eigs)");
  printf("Done.\n");
  return 0;
}
