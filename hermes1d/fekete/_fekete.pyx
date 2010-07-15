cimport cython

from numpy cimport ndarray
from numpy import zeros

def eval_polynomial_orig(coeffs, x):
    r = 0
    n = len(coeffs)
    for i, a in enumerate(coeffs):
        r += a*x**(n-i-1)
    return r

@cython.boundscheck(False)
def eval_polynomial(ndarray[double] coeffs not None, double x):
    cdef double r=0
    cdef unsigned n = len(coeffs)
    cdef unsigned i
    for i in range(n):
        r += coeffs[i]*x**(n-i-1)
    return r

def get_x_phys_orig(x_ref, a, b):
    return (a+b)/2. + x_ref*(b-a)/2.

def get_x_phys(double x_ref, double a, double b):
    return (a+b)/2. + x_ref*(b-a)/2.

def eval_polynomial_array_orig(coeffs, x):
    r = zeros(len(x))
    n = len(coeffs)
    for i, a in enumerate(coeffs):
        r += a*x**(n-i-1)
    return r

@cython.boundscheck(False)
def eval_polynomial_array(ndarray[double] coeffs not None, ndarray[double] x not None):
    cdef unsigned n_coeffs = len(coeffs)
    cdef unsigned n_x = len(x)
    cdef ndarray[double] r=zeros(n_x)
    cdef unsigned i, j
    for j in range(n_x):
        for i in range(n_coeffs):
            r[j] += coeffs[i]*x[j]**(n_coeffs-i-1)
    return r
