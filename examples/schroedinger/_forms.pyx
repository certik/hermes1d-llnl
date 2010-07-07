from hermes1d._hermes1d cimport Mesh
from _hermes_common cimport Matrix

def assemble_schroedinger(Mesh mesh, Matrix A, Matrix B, l=0):
    c_assemble_schroedinger(mesh.thisptr, A.thisptr, B.thisptr, l)
