# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

from hermes_common._hermes_common cimport c_Matrix

cdef extern from "hermes1d.h":

    ctypedef double double4[4]
    ctypedef double double3[3]
    ctypedef int int3[3]
    ctypedef int int2[2]

    cdef cppclass Element:
        double x1, x2
        int p
        int *dof
        double get_solution_value(double x_phys, int comp)
        double get_solution_deriv(double x_phys, int comp)
        void get_coeffs(int sln, int comp, double coeffs[])

    cdef cppclass Mesh:
        Mesh(double a, double b, int n_elem, int p_init, int n_eq, int
                n_sln, int print_banner)
        Mesh(int n_macro_elem, double *pts_array, int *p_array, int *m_array,
             int *div_array, int n_eq, int n_sln, int print_banner)
        void create(double A, double B, int n)
        int get_n_base_elems()
        void set_poly_orders(int poly_order)
        int assign_dofs()
        Element *get_base_elems()
        void set_bc_left_dirichlet(int eq_n, double val)
        void set_bc_right_dirichlet(int eq_n, double val)
        void copy_vector_to_mesh(double *y, int sln)
        void copy_mesh_to_vector(double *y, int sln)
        void plot(char* filename)
        Mesh *replicate()
        void reference_refinement(int start_elem_id, int elem_num)
        int get_n_active_elem()
        int get_n_dof()

    cdef cppclass Linearizer:
        Linearizer(Mesh *mesh)
        void plot_solution(char *out_filename,
                int plotting_elem_subdivision)
        void get_xy_mesh(int comp, int plotting_elem_subdivision,
                double **x, double **y, int *n)

    cdef cppclass Iterator:
        Iterator(Mesh *mesh)
        void reset()
        Element *first_active_element()
        Element *next_active_element()
        Element *last_active_element()

    double calc_error_estimate(int norm, Mesh* mesh, Mesh* mesh_ref,
                       double *err_array, int sln)
    double calc_solution_norm(int norm, Mesh* mesh)

    void adapt(int norm, int adapt_type, double threshold,
               double *err_squared_array,
               Mesh* &mesh, Mesh* &mesh_ref)

    int H1D_L2_ortho_global
    int H1D_H1_ortho_global
    ctypedef void(*ExactFunction)(int n, double x[], double f[], double dfdx[])
    void assemble_projection_matrix_rhs(Mesh *mesh, c_Matrix *A, double *rhs,
            ExactFunction fn, int projection_type) except +
