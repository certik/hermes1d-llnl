from numpy cimport ndarray

def eval_polynomial1(coeffs, x):
    r = 0
    n = len(coeffs)
    for i, a in enumerate(coeffs):
        r += a*x**(n-i-1)
    return r

def eval_polynomial2(ndarray coeffs, double x):
    cdef double r=0
    cdef int n = len(coeffs)
    cdef int i
    for i in range(n):
        r += coeffs[i]*x**(n-i-1)
    return r
