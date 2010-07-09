from hermes1d.h1d_wrapper.h1d_wrapper cimport Mesh
from hermes_common._hermes_common cimport Matrix

def assemble_schroedinger(Mesh mesh, Matrix A, Matrix B, l=0, eqn_type=None):
    cdef int equation_type
    if eqn_type == "R":
        equation_type = c_eqn_type_R
    elif eqn_type == "rR":
        equation_type = c_eqn_type_rR
    else:
        raise ValueError("Unknown equation type")
    c_assemble_schroedinger(mesh.thisptr, A.thisptr, B.thisptr, l,
            equation_type)
