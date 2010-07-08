from hermes1d.h1d_wrapper.h1d_wrapper cimport Mesh
from hermes_common._hermes_common cimport Matrix

def assemble_schroedinger(Mesh mesh, Matrix A, Matrix B, l=0):
    c_assemble_schroedinger(mesh.thisptr, A.thisptr, B.thisptr, l)
