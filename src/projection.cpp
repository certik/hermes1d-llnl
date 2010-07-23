// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "projection.h"
#include "math.h"
#include "discrete.h"

double L2_projection_biform(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for (int i=0; i<num; i++) {
        val += weights[i] * (u[i] * v[i]);
    }
    return val;
}

double H1_projection_biform(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for (int i=0; i<num; i++) {
        val += weights[i] * (u[i] * v[i] + dudx[i] * dvdx[i]);
    }
    return val;
}

// If dx == NULL, don't return the derivative
typedef double(*ExactFunction)(double x, double *dx);

double __f(double x, double *dx)
{
    if (dx != NULL)
        *dx = cos(x);
    return sin(x);
}

ExactFunction _f = __f;

double L2_projection_liform(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for (int i=0; i<num; i++) {
        // Value of the projected function at gauss points:
        double f = _f(x[i], NULL);
        val += weights[i] * (f * v[i]);
    }
    return val;
}

double H1_projection_liform(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for (int i=0; i<num; i++) {
        // Value of the projected function at gauss points:
        double dfdx;
        double f = _f(x[i], &dfdx);
        val += weights[i] * (f * v[i] + dfdx * dvdx[i]);
    }
    return val;
}

void assemble_projection_matrix_rhs(Mesh *mesh, Matrix *A, double *rhs,
        int projection_type)
{
    DiscreteProblem *dp1 = new DiscreteProblem();
    if (projection_type == H1D_L2_ortho_global) {
        dp1->add_matrix_form(0, 0, L2_projection_biform);
        dp1->add_vector_form(0, L2_projection_liform);
    } else if (projection_type == H1D_H1_ortho_global) {
        dp1->add_matrix_form(0, 0, H1_projection_biform);
        dp1->add_vector_form(0, H1_projection_liform);
    } else
        throw std::runtime_error("Unknown projection type");
    int N_dof = mesh->assign_dofs();
    printf("Assembling projection linear system. ndofs: %d\n", N_dof);
    dp1->assemble_matrix_and_vector(mesh, A, rhs);
    printf("  Done assembling.\n");
    delete dp1;
}
