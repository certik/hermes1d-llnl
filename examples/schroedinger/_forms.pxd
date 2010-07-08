from hermes1d.h1d_wrapper.hermes1d cimport Mesh
from _hermes_common cimport c_Matrix

cdef extern from "forms.h":
    void c_assemble_schroedinger "assemble_schroedinger"(Mesh *mesh,
            c_Matrix *A, c_Matrix *B, int l)
