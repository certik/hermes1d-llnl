from hermes1d.h1d_wrapper._hermes1d cimport c_Mesh
from _hermes_common cimport c_Matrix

cdef extern from "forms.h":
    void c_assemble_schroedinger "assemble_schroedinger"(c_Mesh *mesh,
            c_Matrix *A, c_Matrix *B, int l)
