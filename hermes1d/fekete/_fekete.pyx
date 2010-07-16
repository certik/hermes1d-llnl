cimport cython

from numpy cimport ndarray
from numpy import zeros, empty, real, array
from numpy.linalg import solve
from scipy.special.orthogonal import p_roots
from gauss_lobatto_points import points

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

# gauss points+weights on the reference domain
cdef dict _gauss_points_reference = {}

def init_gauss_points():
    print "Precalculating gauss points..."
    for n in range(1, 50):
        x, w = p_roots(n)
        _gauss_points_reference[n] = (real(x), w)
    print "    Done."
init_gauss_points()

@cython.boundscheck(False)
def get_gauss_points(double a, double b, int n):
    cdef double J = (b-a)/2.0
    cdef unsigned n_points, i
    cdef ndarray[double, mode="c"] x, w
    x, w = _gauss_points_reference[n]
    n_points = len(x)
    cdef ndarray[double, mode="c"] x_phys=empty(n_points), w_phys=empty(n_points)
    for i in range(n_points):
        x_phys[i] = J*(x[i]+1) + a
        w_phys[i] = w[i]*J
    return x_phys, w_phys

def eval_poly(n, x, values, a, b):
    """
    Evaluates polynomial at points 'x'.

    The polynomial is defined by 'values' in fekete points, on the interval
    (a, b).
    """
    fekete_points = array(points[len(values)-1])
    coeffs = get_polynomial(fekete_points, values, a, b)
    return eval_polynomial_array(coeffs, x)
