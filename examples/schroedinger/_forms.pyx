from _hermes1d cimport Mesh
from _hermes_common cimport Matrix

def assemble_schroedinger(Mesh mesh, Matrix A, Matrix B):
    c_assemble_schroedinger(mesh.thisptr, A.thisptr, B.thisptr)
