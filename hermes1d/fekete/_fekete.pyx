cimport cython

from numpy cimport ndarray
from numpy import zeros, empty
from numpy.linalg import solve

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

cpdef get_x_phys(double x_ref, double a, double b):
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

def get_polynomial(x, values, a, b):
    """
    Returns the interpolating polynomial's coeffs.

    The len(values) specifies the order and we work in the element <a, b>
    """
    n = len(values)
    A = empty((n, n), dtype="double")
    y = empty((n,), dtype="double")
    assert len(x) == n
    for i in range(n):
        for j in range(n):
            A[i, j] = get_x_phys(x[i], a, b)**(n-j-1)
        y[i] = values[i]
    a = solve(A, y)
    return a
