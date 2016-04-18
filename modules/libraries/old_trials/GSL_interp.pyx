# import numpy as np
# cimport numpy as np
cimport cython
from cython_gsl cimport *

cdef double interpolation_GSL_C(double x, double[::1] data_x, double[::1] data_y):

    cdef int N = data_x.shape[0]

    cdef gsl_interp_accel *acc
    acc = gsl_interp_accel_alloc ()
    cdef gsl_spline *spline
    spline = gsl_spline_alloc (gsl_interp_cspline, N)

    gsl_spline_init (spline, &data_x[0], &data_y[0], N)
    cdef double y = gsl_spline_eval (spline, x, acc)

    gsl_spline_free (spline)
    gsl_interp_accel_free (acc)

    return(y)

def prova_time(double x, double[::1] data_x, double[::1] data_y):
    for i in range(1000000):
        interpolation_GSL_C(x, data_x, data_y)


def interpolation_GSL(double x, double[::1] data_x, double[::1] data_y):
    return(interpolation_GSL_C(x, data_x, data_y))
