cimport cython

from numpy cimport ndarray, double_t

def eval_polynomial1(coeffs, x):
    r = 0
    n = len(coeffs)
    for i, a in enumerate(coeffs):
        r += a*x**(n-i-1)
    return r

@cython.boundscheck(False)
def eval_polynomial2(ndarray[double_t] coeffs not None, double x):
    cdef double r=0
    cdef unsigned n = len(coeffs)
    cdef unsigned i
    for i in range(n):
        r += coeffs[i]*x**(n-i-1)
    return r

def get_x_phys(x_ref, a, b):
    return (a+b)/2. + x_ref*(b-a)/2.
