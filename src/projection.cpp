// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "projection.h"
#include "math.h"

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


double _f(double x)
{
    return sin(x);
}

double _dfdx(double x)
{
    return cos(x);
}

double L2_projection_liform(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for (int i=0; i<num; i++) {
        // Value of the projected function at gauss points:
        double f = _f(x[i]);
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
        double f = _f(x[i]);
        double dfdx = _dfdx(x[i]);
        val += weights[i] * (f * v[i] + dfdx * dvdx[i]);
    }
    return val;
}
