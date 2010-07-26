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

/*
    Returns the values of the function in f[] and derivatives in dfdx[].

    For all physical 'x' in x[]. If dfdx == NULL, don't return the derivative.

    This function is completely general and must work for any number of points
    'n'. Typically, this function is being called on each element, in which
    case the x[] are Gauss integratin points::

        ExactFunction f = <initialize it>;
        double x[10] = <initialize Gauss points>;
        double val[10], dfdx[10];
        f(10, x, val, dfdx);

    If you want to get just a value and derivative at a point, call it like::

        ExactFunction f = <initialize it>;
        double x = 3.15;
        double val, dfdx;
        f(1, &x, &val, &dfdx);

    or if you only need the value::

        ExactFunction f = <initialize it>;
        double x = 3.15;
        double val;
        f(1, &x, &val, NULL);
*/
typedef void(*ExactFunction)(int n, double x[], double f[], double dfdx[]);

void f_sin(int n, double x[], double f[], double dfdx[])
{
    for (int i=0; i<n; i++) {
        f[i] = sin(x[i]);
        if (dfdx != NULL)
            dfdx[i] = cos(x[i]);
    }
}

void f_unimplemented(int n, double x[], double f[], double dfdx[])
{
    error("Internal error: you need to implement _f");
}


ExactFunction _f = f_sin;

double L2_projection_liform(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
    // Value of the projected function at gauss points:
    double f[num];
    _f(num, x, f, NULL);
    double val = 0;
    for (int i=0; i<num; i++) {
        val += weights[i] * (f[i] * v[i]);
    }
    return val;
}

double H1_projection_liform(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
    // Value of the projected function at gauss points:
    double f[num], dfdx[num];
    _f(num, x, f, dfdx);
    double val = 0;
    for (int i=0; i<num; i++) {
        val += weights[i] * (f[i] * v[i] + dfdx[i] * dvdx[i]);
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
