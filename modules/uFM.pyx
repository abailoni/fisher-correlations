
#------------------------
# Importing modules and defining default data:
#------------------------

include "libraries/header.pxi"

include "libraries/CAMB.pxi"

include "libraries/analytical_fcts.pxi"


# It does all the first necessary initialisations:
def init():
    compute_CAMB_spectra()
    compute_survey_DATA()
    store_int1()

#---------------------------------------
# CONVOLVED SPECTRA:
#---------------------------------------

def integral_1(bin1,bin2,name_var,vect_k = np.linspace(0.0001,0.5,10000)):
    N_k = vect_k.shape[0]
    R = vect_k[-1]
    K_samples, P_samples = np.empty(N_k), np.empty(N_k)
    for i in range(N_k):
        K_samples[i] = K(vect_k[i],bin1,bin2)
        if name_var=="spectrum":
            P_samples[i] = zero_spectrum(vect_k[i])
        else:
            P_samples[i] = CAMB_numerical_paramDER(vect_k[i],n_var_import[name_var])

    # Radial FFT:
    return np.sqrt(vol_shell_original(bin1)*vol_shell_original(bin2))/(2*np.pi)**3 *FFTt.radial_convolution(P_samples,K_samples,R)


# Storing and interpolating int1 with GSL:
cdef enum:
    max_N_bins = 50

cdef:
    interpolation_tools integral_1_tools[max_N_bins][max_N_bins][N_vars]
    # interpolation_tools integral_DER[max_N_bins][max_N_bins]

def store_int1():
    print "Computing, storing and interpolating convolved spectra..."
    vect_k = np.linspace(0.0001,0.5,5000)
    import_variables = ["spectrum","h","n_s","Om_b","Om_c"]
    start = time.time()
    for bin in range(N_bins):
        for name_var in import_variables:
            alloc_interp_GSL(vect_k, integral_1(bin,bin,name_var,vect_k), &integral_1_tools[bin][bin][n_var_import[name_var]])
    print "Done! (%g sec.)" %(time.time()-start)

#-----------------------------------------------
# COMPUTING AP term: (simple way for the moment)
#-----------------------------------------------
# k-derivative of the convolved spectrum:
cdef double conv_spectrum_der_k(double k, int bin):
    return gsl_spline_eval_deriv(integral_1_tools[bin][bin][0].spline, k, integral_1_tools[bin][bin][0].acc)

cdef double spectrum_der_k(double k, int bin):
    if flag_tophat==0:
        return zero_spectrum_der_k(k)
    else:
        return conv_spectrum_der_k(k,bin)

#**************************************************************
#**************************************************************
#
# Fisher Matrix computation:
#
#**************************************************************
#**************************************************************

cdef:
    int flag_tophat = 0

#----------------------------------------------
# Reconstruct spectrum and integrals A and B:
#----------------------------------------------
cdef double windowed_zeroSpectrum(double k, int bin):
    return eval_interp_GSL(k, &integral_1_tools[bin][bin][0])

cdef double windowed_numerical_paramDER(double k, int bin, int var): #var [0-3]
    return eval_interp_GSL(k, &integral_1_tools[bin][bin][var+1])


cdef double spectrum(double k, int bin):
    if flag_tophat==0:
        return zero_spectrum(k)
    else:
        return windowed_zeroSpectrum(k,bin)

cdef double numerical_paramDER(double k, int bin, int var): #var [0-3]
    if flag_tophat==0:
        return CAMB_numerical_paramDER(k,var+1)
    else:
        return windowed_numerical_paramDER(k,bin,var)

#--------------------------------------------------------------
# Contructing the final derivatives for the Fisher Matrix:
#--------------------------------------------------------------
# Observed spectrum: (optimized!)
cdef double observed_spectrum(int bin, double k, double mu):
    return Growth_bins[bin]**2*bias_bins[bin]**2 * (1+beta_bins[bin]*mu**2)**2 * spectrum(k,bin)
def observed_spectrum_py(bin,k,mu):
    return observed_spectrum(bin,k,mu)

# Observed terms: (optimized!)
# to avoid division by zero given by windowed_Spectrum with i!=j
cdef double observed_terms(int bin, double k, double mu):
    return Growth_bins[bin]**2 * bias_bins[bin]**2 * (1+beta_bins[bin]*mu**2)**2

def test_growth(bin, k, mu):
    return( Growth_bins[bin]**2 * zero_spectrum(k) )


# h and n_s: (optimized!) (var 0-1)
cdef double der_type_A(int bin,  int var_num, double k,double mu):
    return observed_terms(bin, k, mu) * numerical_paramDER(k,bin,var_num)

# Om_b, Om_c, w_0: (optimized!) (var 2-4)
cdef double der_type_B(int bin, int var_num, double k, double mu):
    cdef double CAMB_term
    if var_num<=3: #Om_b, Om_c
        CAMB_term = numerical_paramDER(k,bin,var_num)
    else:
        CAMB_term = 0
    cdef double beta_term = 2./(1+beta_bins[bin]*mu**2)*( 2*beta_bins[bin]*mu*mu_der(mu,bin,var_num) + mu**2*Beta_der_data[var_num][bin] )
    cdef double k_derivative = spectrum_der_k(k,bin)*k_der(mu, k, bin, var_num)
    return observed_terms(bin,k,mu) * (CAMB_term + k_derivative*CHECK_DER_K) + observed_spectrum(bin, k, mu) * (2*lnG_der_data[var_num][bin] + lnH_der_data[var_num][bin] - 2*lnD_der_data[var_num][bin] + beta_term  )


# # Gamma: (optimized!) (var=6)
# cdef double der_gamma(int bin, double k, double mu):
#     # Pay attention to lnH_der_data that are computed in z_avg....
#     return(observed_spectrum(bin, k, mu) * (2*lnG_der_data[6][bin] + 2./(1+beta_bins[bin]*mu**2)*(mu**2*Beta_der_data[6][bin]) ) )

# Sigma8: (optimized!) (var=5)
cdef double der_sigma8(int bin, double k, double mu):
    return 2*observed_spectrum(bin, k, mu) /ref_values["sigma8"]

# Bias: (9 derivatives) (bad optimized....) (var>=6)
cdef double der_bias(int bin, double k, double mu, int bin_bias) except -1:
    if bin==bin_bias:
        result =observed_spectrum(bin, k, mu) * 2 * (1/bias_bins[bin] - 1./(1+beta_bins[bin]*mu**2)*mu**2 * fnEv(Om_m_z_py,z=z_avg[bin],w_1=ref_values['w_1'],w_0=ref_values['w_0'],Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'])**ref_values['gamma'] /(bias_bins[bin]**2) )
        return(result)
    else:
        return(0.)


# Main routine for derivatives:
cdef double DER(double k, double mu, int var_num, int bin):
    if var_num+1<=2: #h and n_s
        return(der_type_A(bin,var_num,k,mu))
    elif var_num+1<=5: # Om_b, Om_c, w_0:
        return(der_type_B(bin,var_num,k,mu))
    elif var_num+1==6: #sigma8
        return(der_sigma8(bin,k,mu))
    else: #bias
        bin_bias = var_num-N_cosm_vars
        return(der_bias(bin,k,mu,bin_bias))


#--------------------------------------------------------------
#--------------------------------------------------------------
# Computation of the Fisher Matrix: (with GSL)
#--------------------------------------------------------------
#--------------------------------------------------------------

# Integration variables: (GSL)
cdef enum:
    max_alloc_const = 500
cdef:
    size_t MAX_ALLOC = max_alloc_const
    double CHECK_DER_K = 0
    double rel_prec = 1e-8
    double abs_prec = 1e-6
cdef:
    gsl_integration_workspace * W_k
    gsl_integration_workspace * W_mu
W_k = gsl_integration_workspace_alloc(MAX_ALLOC)
W_mu = gsl_integration_workspace_alloc(MAX_ALLOC)

# Integral arguments:
cdef double argument_mu(double mu, void *input): #k, var1, var2, bin
    cdef:
        double *params = <double*> input
        double k = params[0]
        int bin = <int>params[3]
        double eff_term = (n_dens_c[bin]/(n_dens_c[bin]*observed_spectrum(bin,k,mu)+1))**2
    return k*k * DER(k, mu, <int>params[1], bin) * DER(k, mu, <int>params[2], bin) * eff_term

cdef double argument_k(double k, void *input): #var1, var2, bin
    cdef double params_mu[4]
    cdef double *in_params = <double*>input
    cdef gsl_function F_mu
    F_mu.function = &argument_mu
    params_mu[0], params_mu[1], params_mu[2], params_mu[3] = k, in_params[0], in_params[1], in_params[2]
    return eval_integration_GSL(-1., 1., abs_prec, rel_prec, params_mu, W_mu, &F_mu, MAX_ALLOC)


#------------------------
# FISHER MATRIX element:
#------------------------
cdef:
    double k_max_hard = 0.5
    double k_min_hard = 0.001

def fisher_matrix_element(int var1, int var2, int check_AP=0, top_hat=0, double fixed_kmax=0.2):
    global CHECK_DER_K, flag_tophat
    CHECK_DER_K, flag_tophat = check_AP, top_hat
    cdef gsl_function F_k
    F_k.function = &argument_k
    cdef double params[3]
    params[0], params[1] = var1, var2
    cdef double FM = 0
    cdef double[::1] k_max_array = np.empty(N_bins)
    if abs(fixed_kmax)<1e-10: # is zero
        for bin in range(N_bins):
            if (k_max_data[bin]<k_max_hard):
                k_max_array[bin]=k_max_data[bin]
            else:
                k_max_array[bin]=k_max_hard
    else:
        k_max_array = np.array([fixed_kmax]*N_bins)

    # Check if some element of the matrix is zero anyway:
    if var1>N_cosm_vars:
        bin = var1-N_cosm_vars # bin_bias
        params[2]= bin
        FM += vol_shell(bin) /(8*np.pi**2) * eval_integration_GSL(k_min_hard, k_max_array[bin], abs_prec, rel_prec, params, W_k, &F_k, MAX_ALLOC)
    elif var2>N_cosm_vars:
        bin = var2-N_cosm_vars # bin_bias
        params[2]= bin
        FM += vol_shell(bin) /(8*np.pi**2) * eval_integration_GSL(k_min_hard, k_max_array[bin], abs_prec, rel_prec, params, W_k, &F_k, MAX_ALLOC)
    else:
        for bin in range(N_bins):
            params[2]=bin
            FM += vol_shell(bin) /(8*np.pi**2) * eval_integration_GSL(k_min_hard, k_max_array[bin], abs_prec, rel_prec, params, W_k, &F_k, MAX_ALLOC)
    return FM


#------------------------
# Computation of FM:
#------------------------
def FM(int check_AP=0, top_hat=0, FMname="test", double fixed_kmax=0.2):
    FM = np.empty([N_tot_vars,N_tot_vars])
    for var1 in range(N_tot_vars):
        for var2 in range(var1,N_tot_vars):
            if var1>=N_cosm_vars and var2>=N_cosm_vars and var1!=var2:
                FM[var1,var2], FM[var2,var1] = 0., 0.
            else:
                start = time.clock()
                FM[var1,var2]=fisher_matrix_element(var1,var2,check_AP,top_hat,fixed_kmax)
                stop = time.clock()
                FM[var2,var1]=FM[var1,var2]
                print "(%d, %d) --> %g (%g sec.)" %(var1,var2,FM[var1,var2],stop-start)
    np.savetxt("OUTPUT/FMuncorr_%s_AP%d-%dbins-W%d.csv" %(FMname, check_AP, N_bins, top_hat), FM)
    return FM

