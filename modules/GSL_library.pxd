cdef extern from "gsl/gsl_math.h":

  double M_E

  double M_LOG2E

  double M_LOG10E

  double M_SQRT2

  double M_SQRT1_2

  double M_SQRT3

  double M_PI

  double M_PI_2

  double M_PI_4

  double M_SQRTPI

  double M_2_SQRTPI

  double M_1_PI

  double M_2_PI

  double M_LN10

  double M_LN2

  double M_LNPI

  double M_EULER

  int  gsl_isnan(double x) nogil

  int  gsl_isinf(double x) nogil

  int  gsl_finite(double x) nogil

  double  gsl_log1p(double x) nogil

  double  gsl_expm1(double x) nogil

  double  gsl_hypot(double x, double y) nogil

  double  gsl_acosh(double x) nogil

  double  gsl_asinh(double x) nogil

  double  gsl_atanh(double x) nogil

  double  gsl_ldexp(double x, int e) nogil

  double  gsl_frexp(double x, int * e) nogil

  double  gsl_pow_int(double x, int n) nogil

  double  gsl_pow_2(double x) nogil

  double  gsl_pow_3(double x) nogil

  double  gsl_pow_4(double x) nogil

  double  gsl_pow_5(double x) nogil

  double  gsl_pow_6(double x) nogil

  double  gsl_pow_7(double x) nogil

  double  gsl_pow_8(double x) nogil

  double  gsl_pow_9(double x) nogil

  int GSL_SIGN(double x) nogil

  int GSL_IS_ODD(int n) nogil

  int GSL_IS_EVEN(int n) nogil

  double GSL_MAX(double a, double  b) nogil

  double GSL_MIN(double a, double  b) nogil

  double  GSL_MAX_DBL(double a, double b) nogil

  double  GSL_MIN_DBL(double a, double b) nogil

  int  GSL_MAX_INT(int a, int b) nogil

  int  GSL_MIN_INT(int a, int b) nogil

  long double  GSL_MAX_LDBL(long double a, long double b) nogil

  long double  GSL_MIN_LDBL(long double a, long double b) nogil

  int  gsl_fcmp(double x, double y, double epsilon) nogil

  # Definition of an arbitrary function with parameters
  ctypedef struct gsl_function:
    double (* function) (double x, void * params)
    void * params

  double GSL_FN_EVAL(gsl_function * F, double x) nogil

  # Definition of an arbitrary function returning two values, r1, r2
  ctypedef struct gsl_function_fdf:
    double (* f) (double x, void * params) nogil
    double (* df) (double x, void * params) nogil
    void (* fdf) (double x, void * params, double * f, double * df) nogil
    void * params

  double GSL_FN_FDF_EVAL_F(gsl_function_fdf * FDF, double x) nogil

  GSL_FN_FDF_EVAL_DF(gsl_function_fdf * FDF,double x) nogil

  GSL_FN_FDF_EVAL_F_DF(gsl_function_fdf * FDF,double x, double y,double dy) nogil

###########################################
# INTEGRATION:
###########################################

cdef extern from "gsl/gsl_integration.h":

  ctypedef struct gsl_integration_workspace
  ctypedef struct gsl_integration_qaws_table
  ctypedef struct  gsl_integration_qawo_table
  ctypedef struct gsl_integration_cquad_workspace
  cdef enum:
    GSL_INTEG_GAUSS15 = 1
    GSL_INTEG_GAUSS21 = 2
    GSL_INTEG_GAUSS31 = 3
    GSL_INTEG_GAUSS41 = 4
    GSL_INTEG_GAUSS51 = 5
    GSL_INTEG_GAUSS61 = 6
  cdef enum gsl_integration_qawo_enum:
    GSL_INTEG_COSINE, GSL_INTEG_SINE

  gsl_integration_cquad_workspace *  gsl_integration_cquad_workspace_alloc (size_t n) nogil

  void  gsl_integration_cquad_workspace_free (gsl_integration_cquad_workspace * w) nogil

  int gsl_integration_cquad (gsl_function * f, double a, double b, double epsabs, double epsrel, gsl_integration_cquad_workspace * workspace, double * result, double * abserr, size_t * nevals) nogil

  int  gsl_integration_qng(gsl_function *f, double a, double b, double epsabs, double epsrel, double * result, double * abserr, size_t * neval) nogil

  gsl_integration_workspace *  gsl_integration_workspace_alloc(size_t n) nogil

  void  gsl_integration_workspace_free(gsl_integration_workspace * w) nogil

  int  gsl_integration_qag(gsl_function *f, double a, double b, double epsabs, double epsrel, size_t limit, int key, gsl_integration_workspace * workspace, double * result, double * abserr) nogil

  int  gsl_integration_qags(gsl_function * f, double a, double b, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

  int  gsl_integration_qagp(gsl_function * f, double *pts, size_t npts, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

  int  gsl_integration_qagi(gsl_function * f, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

  int  gsl_integration_qagiu(gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

  int  gsl_integration_qagil(gsl_function * f, double b, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

  int  gsl_integration_qawc(gsl_function *f, double a, double b, double c, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr) nogil

  gsl_integration_qaws_table *  gsl_integration_qaws_table_alloc(double alpha, double beta, int mu, int nu) nogil

  int  gsl_integration_qaws_table_set(gsl_integration_qaws_table * t, double alpha, double beta, int mu, int nu) nogil

  void  gsl_integration_qaws_table_free(gsl_integration_qaws_table * t) nogil

  int  gsl_integration_qaws(gsl_function * f, double a, double b, gsl_integration_qaws_table * t, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

  gsl_integration_qawo_table *  gsl_integration_qawo_table_alloc(double omega, double L,  gsl_integration_qawo_enum sine, size_t n) nogil

  int  gsl_integration_qawo_table_set(gsl_integration_qawo_table * t, double omega, double L,  gsl_integration_qawo_enum sine) nogil

  int  gsl_integration_qawo_table_set_length(gsl_integration_qawo_table * t, double L) nogil

  void  gsl_integration_qawo_table_free(gsl_integration_qawo_table * t) nogil

  int  gsl_integration_qawo(gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, gsl_integration_qawo_table * wf, double *result, double *abserr) nogil

  int  gsl_integration_qawf(gsl_function * f, double a, double epsabs, size_t limit, gsl_integration_workspace * workspace, gsl_integration_workspace * cycle_workspace, gsl_integration_qawo_table * wf, double *result, double *abserr) nogil

  double GSL_EMAXITER

  double GSL_EROUND

  double GSL_ESING

  double GSL_EDIVERGE



###########################################
# INTERPOLATION:
###########################################

cdef extern from "gsl/gsl_interp.h":

  ctypedef struct gsl_interp_accel

  ctypedef struct gsl_interp_type
  ctypedef struct gsl_interp

  gsl_interp_type * gsl_interp_linear
  gsl_interp_type * gsl_interp_polynomial
  gsl_interp_type * gsl_interp_cspline
  gsl_interp_type * gsl_interp_cspline_periodic
  gsl_interp_type * gsl_interp_akima
  gsl_interp_type * gsl_interp_akima_periodic

  gsl_interp_accel * gsl_interp_accel_alloc() nogil

  size_t gsl_interp_accel_find(gsl_interp_accel * a,  double x_array[], size_t size, double x) nogil

  int gsl_interp_accel_reset (gsl_interp_accel * a) nogil

  void gsl_interp_accel_free(gsl_interp_accel * a) nogil

  gsl_interp * gsl_interp_alloc( gsl_interp_type * T, size_t n) nogil

  int gsl_interp_init(gsl_interp * obj,  double xa[],  double ya[], size_t size) nogil

  char * gsl_interp_name( gsl_interp * interp) nogil
  unsigned int gsl_interp_min_size( gsl_interp * interp) nogil


  int gsl_interp_eval_e( gsl_interp * obj,
                     double xa[],  double ya[], double x,
                    gsl_interp_accel * a, double * y) nogil

  double gsl_interp_eval( gsl_interp * obj,
                   double xa[],  double ya[], double x,
                  gsl_interp_accel * a) nogil

  int gsl_interp_eval_deriv_e( gsl_interp * obj,
                           double xa[],  double ya[], double x,
                          gsl_interp_accel * a,
                          double * d) nogil

  double gsl_interp_eval_deriv( gsl_interp * obj,
                         double xa[],  double ya[], double x,
                        gsl_interp_accel * a) nogil

  int gsl_interp_eval_deriv2_e( gsl_interp * obj,
                            double xa[],  double ya[], double x,
                           gsl_interp_accel * a,
                           double * d2) nogil

  double gsl_interp_eval_deriv2( gsl_interp * obj,
                          double xa[],  double ya[], double x,
                         gsl_interp_accel * a) nogil

  int gsl_interp_eval_integ_e( gsl_interp * obj,
                           double xa[],  double ya[],
                          double a, double b,
                          gsl_interp_accel * acc,
                          double * result) nogil

  double gsl_interp_eval_integ( gsl_interp * obj,
                         double xa[],  double ya[],
                        double a, double b,
                        gsl_interp_accel * acc) nogil

  void gsl_interp_free(gsl_interp * interp) nogil

  size_t gsl_interp_bsearch( double x_array[], double x,
                            size_t index_lo, size_t index_hi) nogil



cdef extern from "gsl/gsl_spline.h":
  ctypedef struct gsl_spline

  gsl_spline * gsl_spline_alloc( gsl_interp_type * T, size_t size) nogil

  int gsl_spline_init(gsl_spline * spline,  double xa[],  double ya[], size_t size) nogil


  int gsl_spline_eval_e( gsl_spline * spline, double x,
                    gsl_interp_accel * a, double * y) nogil

  double gsl_spline_eval( gsl_spline * spline, double x, gsl_interp_accel * a) nogil

  int gsl_spline_eval_deriv_e( gsl_spline * spline, double x,
                          gsl_interp_accel * a, double * y) nogil

  double gsl_spline_eval_deriv( gsl_spline * spline, double x, gsl_interp_accel * a) nogil

  int gsl_spline_eval_deriv2_e( gsl_spline * spline, double x,
                           gsl_interp_accel * a, double * y) nogil

  double gsl_spline_eval_deriv2( gsl_spline * spline, double x,
                         gsl_interp_accel * a) nogil

  int gsl_spline_eval_integ_e( gsl_spline * spline, double a, double b,
                          gsl_interp_accel * acc, double * y) nogil

  double gsl_spline_eval_integ( gsl_spline * spline, double a, double b,
                        gsl_interp_accel * acc) nogil

  void gsl_spline_free(gsl_spline * spline) nogil


cdef extern from "gsl/gsl_spline2d.h":
  ctypedef struct gsl_spline2d
  ctypedef struct gsl_interp2d_type

  gsl_interp2d_type * gsl_interp2d_bilinear
  gsl_interp2d_type * gsl_interp2d_bicubic

  gsl_spline2d * gsl_spline2d_alloc (gsl_interp2d_type * T, size_t xsize, size_t ysize) nogil

  int gsl_spline2d_init (gsl_spline2d * spline, double xa[], double ya[], double za[], size_t xsize, size_t ysize) nogil

  double gsl_spline2d_eval (gsl_spline2d * spline, double x, double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc) nogil
  int gsl_spline2d_eval_e (gsl_spline2d * spline, double x, double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * z) nogil

  double gsl_spline2d_eval_deriv_x (gsl_spline2d * spline, double x, double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc) nogil
  int gsl_spline2d_eval_deriv_x_e (gsl_spline2d * spline, double x, double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d) nogil

  double gsl_spline2d_eval_deriv_y (gsl_spline2d * spline, double x, double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc) nogil
  int gsl_spline2d_eval_deriv_y_e (gsl_spline2d * spline, double x, double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d) nogil

  void gsl_spline2d_free (gsl_spline2d * spline) nogil



cdef extern from "gsl/gsl_errno.h":
  ctypedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno)
  gsl_error_handler_t* gsl_set_error_handler_off()


