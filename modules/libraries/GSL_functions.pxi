#**************************************************************
#**************************************************************
#
# GSL routines written in Cython:
#
#**************************************************************
#**************************************************************

############################################
# GSL INTERPOLATION:
############################################
cdef struct interpolation_tools:
    gsl_spline* spline
    gsl_interp_accel* acc

@cython.boundscheck(False) # turn off array bounds check
@cython.wraparound(False) # turn off negative indices ([-1,-1])
cdef void alloc_interp_GSL(double[::1] data_x, double[::1] data_y, interpolation_tools *interp_data):
    cdef int N = data_x.shape[0]
    interp_data.acc = gsl_interp_accel_alloc()
    interp_data.spline = gsl_spline_alloc(gsl_interp_cspline, N)
    # Computation spline coeff.:
    gsl_spline_init(interp_data.spline, &data_x[0], &data_y[0], N)
    return

cdef double eval_interp_GSL(double x, interpolation_tools *interp_data):
    return(gsl_spline_eval(interp_data.spline, x, interp_data.acc))

cdef void free_interp_GSL(interpolation_tools *interp_data):
    gsl_spline_free (interp_data.spline)
    gsl_interp_accel_free (interp_data.acc)
    return

# This routine does everything at once:
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double interpolation_GSL(double x, double[::1] data_x, double[::1] data_y):
    cdef int N = data_x.shape[0]
    # Allocating:
    cdef gsl_interp_accel *acc
    acc = gsl_interp_accel_alloc ()
    cdef gsl_spline *spline
    spline = gsl_spline_alloc (gsl_interp_cspline, N)
    # Computation:
    gsl_spline_init (spline, &data_x[0], &data_y[0], N)
    cdef double y = gsl_spline_eval (spline, x, acc)
    gsl_spline_free (spline)
    gsl_interp_accel_free (acc)
    return(y)


###############################################
# GSL INTEGRATION:
###############################################
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double eval_integration_GSL_noadp(double a, double b, double abs_prec, double rel_prec, void *params, gsl_function *F,  size_t *n_eval):
    F.params = params
    cdef double C_results[2]
    #cdef double[::1] results = C_results
    gsl_integration_qng(F,a,b, abs_prec, rel_prec, &C_results[0], &C_results[1], n_eval)
    #gsl_integration_qag(F, a, b, abs_prec, rel_prec, max_alloc, GSL_INTEG_GAUSS15, W, &C_results[0], &C_results[1])
    return(C_results[0])

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double eval_integration_GSL(double a, double b, double abs_prec, double rel_prec, void *params, gsl_integration_workspace *W, gsl_function *F, size_t max_alloc):
    F.params = params
    cdef double C_results[2]
    #cdef double[::1] results = C_results
    cdef int status = gsl_integration_qag(F, a, b, abs_prec, rel_prec, max_alloc, GSL_INTEG_GAUSS15, W, &C_results[0], &C_results[1])
    if status!=0:
        print "#!#",
    return(C_results[0])

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double eval_integration_GSL_61(double a, double b, double abs_prec, double rel_prec, void *params, gsl_integration_workspace *W, gsl_function *F, size_t max_alloc):
    F.params = params
    cdef double C_results[2]
    #cdef double[::1] results = C_results
    gsl_integration_qag(F, a, b, abs_prec, rel_prec, max_alloc, GSL_INTEG_GAUSS61, W, &C_results[0], &C_results[1])
    return(C_results[0])

# n_points includes extremes
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double eval_integration_GSL_sing(double a, double b, double *singularities, size_t n_points, double abs_prec, double rel_prec, void *params, gsl_integration_workspace *W, gsl_function *F, size_t max_alloc):
    F.params = params
    cdef:
        double points[100] #!!!!
        double C_results[2]
        #double[::1] results = C_results
    points[0], points[n_points-1] = a, b
    for i in range(1,n_points-1):
        points[i]=singularities[i-1]
    gsl_integration_qagp(F, &points[0], n_points, abs_prec, rel_prec, max_alloc, W, &C_results[0], &C_results[1])
    return(C_results[0])


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double eval_integ_GSL_cquad(double a, double b, double abs_prec, double rel_prec, void *params, gsl_integration_cquad_workspace *W, gsl_function *F, size_t *n_eval):
    F.params = params
    cdef double C_results[2]
    cdef double[::1] results = C_results
    gsl_integration_cquad(F, a, b, abs_prec, rel_prec, W, &C_results[0], &C_results[1], n_eval)
    return(C_results[0])


############################################
# Manual integration (trapezi algorithm)
############################################

'''
This variant starts already with a minimum dx
INPUT int_trapezi():
 - array con estremi
 - precisione relativa per modalità 1 // passo dx per modalità 2
 - eventuali parametri funzione [puntatore void]
 - funzione
 - number of starting intervals
 - modalità: 1 o 2 [int] --> no longer implemented
'''

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double int_trapezi(double x_min, double x_max, double precisione, void *params, double (*fun) (double,void*), int starting_intervals):
    # CASO IN CUI SI VUOLE UNA CERTA PRECISIONE:
    # Calcolo f(x_max) e f(x_min) in modo da dare la prima stima dell'integrale semplicemente come area di un trapezio di altezza x_max-x_min:
    cdef double sum = (fun(x_max,params)+fun(x_min,params))/2.
    cdef double new_sum, total_sum, h=(x_max-x_min)/starting_intervals
    # Compute the first sum:
    for i in range(starting_intervals+1):
        sum=sum+fun(x_min+i*h,params)
    sum=sum*h
    h=h/2
    # Continuo finché non raggiungo la precisione:
    cdef:
        double old_sum=sum
        int intervals=starting_intervals
    for contatore in range(1000):
        # Calcolo i nuovi termini della sommatoria:
        new_sum = 0
        if contatore>5:
            print total_sum
        for i in range(intervals):
            new_sum = new_sum + fun(x_min+h+2*h*i,params)
        new_sum = new_sum*h
        # Sommo i nuovi contributi a quelli vecchi (dato che h è dimezzato quelli vecchi dovranno essere divisi per 2):
        total_sum = new_sum+old_sum/2
        h = h/2
        intervals = intervals*2
        if ((abs(total_sum-old_sum)/total_sum)<=precisione):
            return(total_sum)
        old_sum=total_sum
    print "ERROR trapezi!"
    return(-999999)
    #//CASO IN CUI h E' BEN DEFINITO:
    #} else {
    #    double sum=(fun(x_max,params)+fun(x_min,params))/2;
    #    double h=precisione;
    #    for(double x=estremi[0]+h;x<estremi[1];x+=(h)) {
    #        sum+=fun(x,params);
    #    }
    #    return(sum*=h);
    #}
