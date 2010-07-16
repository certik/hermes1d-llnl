cimport cython

from numpy cimport ndarray
from numpy import zeros, empty, real
from numpy.linalg import solve
from scipy.special.orthogonal import p_roots

def eval_polynomial_orig(coeffs, x):
    r = 0
    n = len(coeffs)
    for i, a in enumerate(coeffs):
        r += a*x**(n-i-1)
    return r

@cython.boundscheck(False)
def eval_polynomial(ndarray[double, mode="c"] coeffs not None, double x):
    cdef double r=0
    cdef unsigned n = len(coeffs)
    cdef unsigned i
    for i in range(n):
        r += coeffs[i]*x**(n-i-1)
    return r

def get_x_phys_orig(x_ref, a, b):
    return (a+b)/2. + x_ref*(b-a)/2.

cpdef double get_x_phys(double x_ref, double a, double b):
    return (a+b)/2. + x_ref*(b-a)/2.

def eval_polynomial_array_orig(coeffs, x):
    r = zeros(len(x))
    n = len(coeffs)
    for i, a in enumerate(coeffs):
        r += a*x**(n-i-1)
    return r

@cython.boundscheck(False)
def eval_polynomial_array(ndarray[double, mode="c"] coeffs not None, ndarray[double, mode="c"] x not None):
    cdef unsigned n_coeffs = len(coeffs)
    cdef unsigned n_x = len(x)
    cdef ndarray[double, mode="c"] r=zeros(n_x)
    cdef unsigned i, j
    for j in range(n_x):
        for i in range(n_coeffs):
            r[j] += coeffs[i]*x[j]**(n_coeffs-i-1)
    return r

def get_polynomial_orig(x, values, a, b):
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

@cython.boundscheck(False)
def get_polynomial(ndarray[double, mode="c"] x not None,
        ndarray[double, mode="c"] values not None,
        double a, double b):
    """
    Returns the interpolating polynomial's coeffs.

    The len(values) specifies the order and we work in the element <a, b>
    """
    cdef unsigned n = len(values)
    cdef ndarray[double, ndim=2, mode="c"] A = empty((n, n), dtype="double")
    cdef ndarray[double, mode="c"] y = empty((n,), dtype="double")
    assert len(x) == n
    cdef unsigned i, j
    for i in range(n):
        for j in range(n):
            A[i, j] = get_x_phys(x[i], a, b)**(n-j-1)
        y[i] = values[i]
    r = solve(A, y)
    return r

@cython.boundscheck(False)
def int_f2(ndarray[double, mode="c"] w not None,
        ndarray[double, mode="c"] values not None):
    cdef double r=0
    cdef unsigned n = len(w)
    cdef unsigned i
    for i in range(n):
        r += (values[i]**2)*w[i]
    return r

__p_roots = {}
def _p_roots(n):
    if n not in __p_roots:
        x, w = p_roots(n)
        __p_roots[n] = (real(x), w)
    return __p_roots[n]

def get_gauss_points(double a, double b, int n):
    cdef double J = (b-a)/2.0
    cdef ndarray[double] x, w
    x, w = _p_roots(n)
    x_phys = J*(x+1) + a
    w = w*J
    return x_phys, w
