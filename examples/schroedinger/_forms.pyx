from hermes1d.h1d_wrapper.h1d_wrapper cimport Mesh
from hermes_common._hermes_common cimport Matrix

from numpy import array, empty, linspace
from scipy.interpolate import InterpolatedUnivariateSpline
from h5py import File

s = None

def assemble_schroedinger(Mesh mesh, Matrix A, Matrix B, l=0, Z=1,
        eqn_type=None):
    cdef int equation_type
    if eqn_type == "R":
        equation_type = c_eqn_type_R
    elif eqn_type == "rR":
        equation_type = c_eqn_type_rR
    else:
        raise ValueError("Unknown equation type")

    f = File("data.hdf5")
    r = array(f["/dft/r"])
    zeff = array(f["/dft/zeff"])
    global s
    s = InterpolatedUnivariateSpline(r, zeff)

    c_assemble_schroedinger(mesh.thisptr, A.thisptr, B.thisptr, l, Z,
            equation_type)

cdef api double Z_eff(double x):
    return s(x)
