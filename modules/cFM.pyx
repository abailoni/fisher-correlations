###########
# Branch:
###########
# - modified Trace that set k_max for the lower index (commented there is also the higher case)
# - FFT used





# Importing modules and defining default data:

include "libraries/header.pxi"

# include "libraries/import_CLASS.pxi"

include "libraries/CAMB.pxi"
include "libraries/analytical_fcts.pxi"
include "libraries/import_AP_integrals.pxi"

# It does all the first necessary things:
def init():
    compute_CAMB_spectra()
    # import_zero_spectrum_der_k()
    compute_survey_DATA()
    store_conv_spectra()
    # import_AP_int()
    # store_int1()



#---------------------------------------
# COMPUTING INTEGRAL 1:
#---------------------------------------

# Compute at once all the correlations using numpy:
def FFTconvoloved_spectra(name_var,vect_k):
    N_k = vect_k.shape[0]
    R = vect_k[-1]
    K_samples, P_samples, Kder_samples, Pder_samples = np.empty(N_k), np.empty(N_k), np.empty(N_k), np.empty(N_k)

    # Construct all bin pairs: (not so pythonic...)
    count=0
    bin_pairs = np.zeros((2, N_bins*(N_bins+1)/2), dtype=np.int)
    for bin1 in range(N_bins):
        for bin2 in range(bin1,N_bins):
            bin_pairs[0,count], bin_pairs[1,count] = bin1, bin2
            count+=1

    # Compute K, P, derP:
    K_samples = K_FFTconvolution(bin_pairs,vect_k)
    if name_var!="spectrum":
        # The best would be to use the num. derivative
        Pder_samples = np.tile(CAMB_numerical_paramDER_py(vect_k, n_var_import[name_var]) , (bin_pairs.shape[1],1))
    else:
        P_samples = np.tile(zero_spectrum_py(vect_k), (bin_pairs.shape[1],1))
        # Kder_samples = (K_FFTconvolution(bin_pairs,vect_k,**{name_var: ref_values[name_var]+epsilon}) - K_FFTconvolution(bin_pairs,vect_k,**{name_var: ref_values[name_var]-epsilon}) ) / (2*epsilon)

    # Radial convolutions:
    if name_var=="spectrum":
        return FFTt.radial_convolution(P_samples,K_samples,R)
    else:
        return FFTt.radial_convolution(Pder_samples,K_samples,R)

def store_conv_spectra():
    print "\nComputing, storing and interpolating convolved spectra..."
    vect_k = np.linspace(0.0,1.4,7000)
    import_variables = ["spectrum","h","n_s","Om_b","Om_c"]

    # Convolve:
    start = time.clock()
    fft1 = {}
    for name_var in import_variables:
        fft1[name_var] = FFTconvoloved_spectra(name_var, vect_k)

    # Interpolate:
    count=0
    for bin1 in range(N_bins):
        for bin2 in range(bin1,N_bins):
            for name_var in import_variables:
                alloc_interp_GSL(vect_k, fft1[name_var][count,:], &integral_1_tools[bin1][bin2][n_var_import[name_var]])
                # if name_var!="spectrum":
                #     alloc_interp_GSL(vect_k, fft2[name_var][count,:], &integral_2_NEW_tools[bin1][bin2][n_var_import[name_var]])
            count+=1
    print "--> Done! (%g sec.)\n" %(time.clock()-start)


def integral_1(bin1,bin2,name_var,vect_k = np.linspace(0.0,0.5,10000)):
    N_k = vect_k.shape[0]
    R = vect_k[-1]
    K_samples, P_samples = np.empty(N_k), np.empty(N_k)

    # First values at k=0:
    K_samples[0] = 1.
    P_samples[0] = 0.

    # Next values:
    for i in range(1,N_k):
        K_samples[i] = K(vect_k[i],bin1,bin2)
        if name_var=="spectrum":
            P_samples[i] = zero_spectrum(vect_k[i])
        else:
            P_samples[i] = CAMB_numerical_paramDER(vect_k[i],n_var_import[name_var])
    # Radial FFT:
    return np.sqrt(vol_shell_original(bin1)*vol_shell_original(bin2))/(2*np.pi)**3 *FFTt.radial_convolution(P_samples,K_samples,R)


def test_integral_1(bin1,bin2,name_var,vect_k = np.linspace(0.0,0.5,10000)):
    N_k = vect_k.shape[0]
    R = vect_k[-1]
    K_samples, P_samples = np.empty(N_k), np.empty(N_k)
    for i in range(N_k):
        K_samples[i] = K(vect_k[i],bin1,bin2)
        if name_var=="spectrum":
            if vect_k[i]>=0.05 and vect_k[i]<0.2:
                P_samples[i] = 1e4
            else:
                P_samples[i] = 0.
    # Radial FFT:
    return np.sqrt(vol_shell_original(bin1)*vol_shell_original(bin2))/(2*np.pi)**3 *FFTt.radial_convolution(P_samples,K_samples,R)


# Storing and interpolating int1 with GSL:
cdef enum:
    max_N_bins = 250

cdef:
    interpolation_tools integral_1_tools[max_N_bins][max_N_bins][N_vars]
    interpolation_tools integral_1_tools_OLD[max_N_bins][max_N_bins][N_vars]
    # interpolation_tools integral_2_NEW_tools[max_N_bins][max_N_bins][N_vars]
    # interpolation_tools integral_DER[max_N_bins][max_N_bins]



def store_int1():
    print "\nComputing, storing and interpolating convolved spectra..."
    vect_k = np.linspace(0.0,1.4,7000)
    import_variables = ["spectrum","h","n_s","Om_b","Om_c"]
    start = time.clock()
    for bin1 in range(N_bins):
        for bin2 in range(bin1,N_bins):
            for name_var in import_variables:
                alloc_interp_GSL(vect_k, integral_1(bin1,bin2,name_var,vect_k), &integral_1_tools_OLD[bin1][bin2][n_var_import[name_var]])
    print "--> Done! (%g sec.)\n" %(time.clock()-start)

def store_int1_test():
    print "Computing, storing and interpolating convolved spectra..."
    vect_k = np.linspace(0.,1.4,7000)
    # vect_k = np.arange(0.,0.401,1e-4)
    import_variables = ["spectrum","h","n_s","Om_b","Om_c"]
    start = time.time()
    for bin1 in range(N_bins):
        for bin2 in range(bin1,N_bins):
            for name_var in import_variables:
                alloc_interp_GSL(vect_k, test_integral_1(bin1,bin2,name_var,vect_k), &integral_1_tools[bin1][bin2][n_var_import[name_var]])
    print "Done! (%g sec.)" %(time.time()-start)

#-----------------------------------------------
# COMPUTING AP term: (simple way for the moment) ---> the only way
#-----------------------------------------------

#  WRONG! No longer used...
# k-derivative of the convolved spectrum:
cdef double conv_spectrum_der_k(double k, int bin1, int bin2):
    return gsl_spline_eval_deriv(integral_1_tools[bin1][bin2][0].spline, k, integral_1_tools[bin1][bin2][0].acc)

cdef double spectrum_der_k(double k, int bin1, int bin2):
    if "windowFun" in typeFM:
        return conv_spectrum_der_k(k,bin1,bin2) # WRONG
    else:
        return zero_spectrum_der_k(k)

def spectrum_der_k_py(k,bin1,bin2):
    return spectrum_der_k(k,bin1,bin2)

#**************************************************************
#**************************************************************
#
# Fisher Matrix computation:
#
#**************************************************************
#**************************************************************

#----------------------------------------------
# Reconstruct spectrum and integrals A and B:
#----------------------------------------------
cdef double windowed_zeroSpectrum(double k, int bin1, int bin2):
    if k<k_min:
        k=k_min
    if bin1<=bin2: #should not be necessary, but...
        return eval_interp_GSL(k, &integral_1_tools[bin1][bin2][0])
    else:
        return eval_interp_GSL(k, &integral_1_tools[bin2][bin1][0])
def windowed_zeroSpectrum_py(k,bin1,bin2):
    """ Vectorized """
    if numpy_check(k):
        results = np.empty(k.shape)
        for i, kval in enumerate(k):
            results[i] = windowed_zeroSpectrum(kval,bin1,bin2)
        return results
    else:
        return windowed_zeroSpectrum(k,bin1,bin2)

# def windowed_zeroSpectrum_NEW(double k, int bin1, int bin2):
#     if k<k_min:
#         k=k_min
#     if bin1<=bin2: #should not be necessary, but...
#         return eval_interp_GSL(k, &integral_1_NEW_tools[bin1][bin2][0])
#     else:
#         return eval_interp_GSL(k, &integral_1_NEW_tools[bin2][bin1][0])

# def der_param_R_vol(double k, int bin1, int bin2, int var): #var [0-3]
#     if k<k_min:
#         k=k_min
#     if bin1<=bin2: #should not be necessary, but...
#         return eval_interp_GSL(k, &integral_2_NEW_tools[bin1][bin2][var+1])
#     else:
#         return eval_interp_GSL(k, &integral_2_NEW_tools[bin2][bin1][var1])

cdef double windowed_numerical_paramDER(double k, int bin1, int bin2, int var): #var [0-3]
    if k<k_min:
        k=k_min
    if bin1<=bin2:
        return eval_interp_GSL(k, &integral_1_tools[bin1][bin2][var+1])
    else:
        return eval_interp_GSL(k, &integral_1_tools[bin2][bin1][var+1])

def windowed_numerical_paramDER_py(k,bin1,bin2,var):
    """ Vectorized """
    if numpy_check(k):
        results = np.empty(k.shape)
        for i, kval in enumerate(k):
            results[i] = windowed_numerical_paramDER(kval,bin1,bin2,var)
        return results
    else:
        return windowed_numerical_paramDER(k,bin1,bin2,var)

cdef double spectrum(double k, int bin1, int bin2):
    if "correlations" not in typeFM and bin1!=bin2:
        return 0.
    if "windowFun" in typeFM:
        return windowed_zeroSpectrum(k,bin1,bin2)
    else:
        return zero_spectrum(k)

# THIS REPETITION OF CODE IS AWFUL....
def spectrum_py(k,bin1,bin2):
    """ Vectorized """
    if "correlations" not in typeFM and bin1!=bin2:
        if numpy_check(k):
            return np.zeros(k.shape)
        else:
            return 0.
    if "windowFun" in typeFM:
        return windowed_zeroSpectrum_py(k,bin1,bin2)
    else:
        return zero_spectrum_py(k)



cdef double numerical_paramDER(double k, int bin1, int bin2, int var): #var [0-3]
    if "correlations" not in typeFM and bin1!=bin2:
        return 0.
    if "windowFun" in typeFM:
        return windowed_numerical_paramDER(k,bin1,bin2,var)
    else:
        return CAMB_numerical_paramDER(k,var+1)

def numerical_paramDER_py(k,bin1,bin2,var):
    """ Vectorized """
    if "correlations" not in typeFM and bin1!=bin2:
        if numpy_check(k):
            return np.zeros(k.shape)
        else:
            return 0.
    if "windowFun" in typeFM:
        return windowed_numerical_paramDER_py(k,bin1,bin2,var)
    else:
        return CAMB_numerical_paramDER_py(k,var+1)
#--------------------------------------------------------------
# Contructing the final derivatives for the Fisher Matrix:
#--------------------------------------------------------------
# Redshift distortion: (argument exponential)
cdef double redshift_exp(int bin1, int bin2, double k, double mu):
    return exp( -0.5*gsl_pow_2(k*mu*redshift_sigma) * ( gsl_pow_2((1+z_avg[bin1])/H_bins[bin1]) + gsl_pow_2((1+z_avg[bin2])/H_bins[bin2]) ) )


# Observed spectrum: (optimized!)
cdef double observed_spectrum(int bin1, int bin2, double k, double mu):
    return Growth_bins[bin1]*Growth_bins[bin2] * bias_bins[bin1]*bias_bins[bin2]* (1+beta_bins[bin1]*mu**2)*(1+beta_bins[bin2]*mu**2) * spectrum(k,bin1,bin2) * redshift_exp(bin1,bin2,k,mu)

# NOT VECTORIZED...
def observed_spectrum_py(bin1,bin2,k,mu):
    """ Vectorized """
    return Growth_bins[bin1]*Growth_bins[bin2] * bias_bins[bin1]*bias_bins[bin2]* (1+beta_bins[bin1]*mu**2)*(1+beta_bins[bin2]*mu**2) * spectrum_py(k,bin1,bin2)* redshift_exp(bin1,bin2,k,mu)

def growth_spectrum_py(bin1,bin2,k,typeFM_input):
    """ Vectorized """
    global typeFM
    typeFM = typeFM_input
    return Growth_bins[bin1]*Growth_bins[bin2] * spectrum_py(k,bin1,bin2)

# Observed terms: (optimized!)
# to avoid division by zero given by windowed_Spectrum with i!=j
cdef double observed_terms(int bin1, int bin2, double k, double mu):
    if "correlations" not in typeFM and bin1!=bin2:
        return 0.
    else:
        return Growth_bins[bin1]*Growth_bins[bin2] * bias_bins[bin1]*bias_bins[bin2]* (1+beta_bins[bin1]*mu**2)*(1+beta_bins[bin2]*mu**2) * redshift_exp(bin1,bin2,k,mu)

# h and n_s: (optimized!) (var 0-1)
cdef double der_type_A(int bin1, int bin2, double k, double mu, int var_num):
    return  observed_terms(bin1, bin2, k, mu)*numerical_paramDER(k, bin1, bin2, var_num)




# Om_b, Om_c, w_0: (optimized!) (var 2-4)
cdef double der_type_B(int bin1, int bin2, double k, double mu, int var_num):
    cdef double CLASS_term
    if var_num<=3: #Om_b, Om_c
        CLASS_term = numerical_paramDER(k, bin1, bin2, var_num)
    else:
        CLASS_term = 0

    check_AP = 1 if AP_flag else 0

    cdef double beta_term = 1./(1+beta_bins[bin1]*mu**2)*( check_AP * 2*beta_bins[bin1]*mu*mu_der(mu,bin1,var_num) + mu**2*Beta_der_data[var_num][bin1] ) +  1./(1+beta_bins[bin2]*mu**2) * ( check_AP * 2*beta_bins[bin2]*mu*mu_der(mu,bin2,var_num)+ mu**2*Beta_der_data[var_num][bin2] )
    # cdef double beta_term = mu_num_term(k, mu, bin1, var_num)


    # Redshift distortion term:
    redshift_distortion = k*mu * gsl_pow_2(redshift_sigma*(1+z_avg[bin1])/H_bins[bin1]) * ( k*mu*lnH_der_data[var_num][bin1] -check_AP*mu*k_der(mu,k,bin1,var_num) -check_AP*k*mu_der(mu,bin1,var_num) + k*mu*lnH_der_data[var_num][bin2] -check_AP*mu*k_der(mu,k,bin2,var_num) -check_AP*k*mu_der(mu,bin2,var_num))


    # AP TERM for k:
    # check_AP = 0. # PUT ALWAYS TO ZERO
    # cdef double AP_term = check_AP * Pk_AP_num_der(k, mu, bin1, var_num)

    # OLD  with dP/dk:
    # cdef double AP_term = check_AP * spectrum_der_k(k,bin1,bin2) * sqrt(k_der(mu,k,bin1,var_num)*k_der(mu,k,bin2,var_num))

    cdef double AP_k_term = 0.
    if AP_flag:
        if "windowFun" not in typeFM:
            AP_k_term = observed_terms(bin1, bin2, k, mu) * check_AP * spectrum_der_k(k,bin1,bin2) * k_der(mu,k,bin1,var_num)
        else:
            var = names_vars[var_num]
            if bin1==bin2:
                AP_k_term = observed_terms(bin1, bin2, k, mu) * check_AP * AP_integral[var][bin1](k,mu)
            else:
                AP_k_term = observed_terms(bin1, bin2, k, mu) * check_AP * AP_integral_corr[var][bin1](k,mu) + observed_spectrum(bin1, bin2, k, mu) * 1/2.*(lnH_der_data[var_num][bin1] - 2*lnD_der_data[var_num][bin1] - lnH_der_data[var_num][bin2] + 2*lnD_der_data[var_num][bin2])
    # cdef double AP_k_term = check_AP * spectrum_der_k(k,bin1,bin2) * k_der(mu,k,bin1,var_num) if bin1==bin2 else 0.

    return observed_terms(bin1, bin2, k, mu)*CLASS_term + AP_k_term + observed_spectrum(bin1, bin2, k, mu) * (lnG_der_data[var_num][bin1]+lnG_der_data[var_num][bin2] + 1/2.*(lnH_der_data[var_num][bin1] - 2*lnD_der_data[var_num][bin1]) + 1/2.*(lnH_der_data[var_num][bin2] - 2*lnD_der_data[var_num][bin2]) + beta_term + zDist_flag*redshift_distortion)



# # Gamma: (optimized!) (var=6)
# cdef double der_gamma(int bin1, int bin2, double k, double mu):
#     # Pay attention to lnH_der_data that are computed in z_avg....
#     return(observed_spectrum(bin1, bin2, k, mu) * (lnG_der_data[6][bin1]+lnG_der_data[6][bin2] + 1./(1+beta_bins[bin1]*mu**2)*(mu**2*Beta_der_data[6][bin1]) + 1./(1+beta_bins[bin2]*mu**2)*(mu**2*Beta_der_data[6][bin2])) )

# def ciao():
#     k1 = 0.1
#     mu1 = 0.5
#     mu = 0.5
#     integral_phi1 = Lambda_Ev(phi1_intregral_data, k1, mu1, mu, lnH_der_data[num_var][bin1], lnD_der_data[num_var][bin1])
#     integral_

# Sigma8: (optimized!) (var=5)
cdef double der_sigma8(int bin1, int bin2, double k, double mu):
    return 2*observed_spectrum(bin1, bin2, k, mu) /ref_values["sigma8"]


cdef double bias_term(double mu, int bin):
    return 1/bias_bins[bin] - 1./(1+beta_bins[bin]*mu**2)*mu**2 * unfun_Ev(Om_m_z_data,z_avg[bin])**ref_values['gamma'] /(bias_bins[bin]**2)


# Bias: (9 derivatives) (bad optimized....) (var>=6)
cdef double der_bias(int bin1, int bin2, double k, double mu, int bin_bias):
    if bin1==bin2 and bin1==bin_bias:
        return observed_spectrum(bin1, bin2, k, mu) * 2*bias_term(mu,bin_bias)
    elif bin1==bin_bias or bin2==bin_bias:
        return observed_spectrum(bin1, bin2, k, mu) * bias_term(mu,bin_bias)
    else:
        return 0.



#--------------------------------------------------------------
# Constructing matrices and find the final Trace: (numpy+cython)
# (optimized with memory-views)
#--------------------------------------------------------------

# Compute inverse covariance matrix:
def inverse_matrix_C(double k, double mu):
    cdef np.intp_t bin1, bin2
    for bin1 in range(N_bins):
        stop_bin2  = bin1+CORR_BINS+1 if bin1+CORR_BINS+1<N_bins else N_bins
        for bin2 in range(bin1,stop_bin2):
            C_v[bin1,bin2]= observed_spectrum(bin1, bin2, k, mu) + N_v[bin1,bin2]
            C_v[bin2,bin1]=C_v[bin1,bin2]
    return(np.linalg.inv(C))

# Main routine for derivatives:
cdef double DER(double k, double mu, int var_num, int bin1, int bin2):
    if var_num+1<=2: #h and n_s
        return der_type_A(bin1,bin2,k,mu,var_num)
    elif var_num+1<=5: # Om_b, Om_c, w_0:
        return der_type_B(bin1,bin2,k,mu,var_num)
    elif var_num+1==6: #sigma8
        return der_sigma8(bin1,bin2,k,mu)
    else: #bias
        bin_bias = var_num-N_cosm_vars_max
        return der_bias(bin1,bin2,k,mu,bin_bias)


# Compute matrix of derivatives: (var_num from 0 to 7+N_bins-1)
cdef void derivative_matrices(double k, double mu, int var_num, double[:,::1] P_der_matrix):
    # P_der_matrix = np.zeros((N_bins,N_bins))
    cdef np.intp_t bin1, bin2
    for bin1 in range(N_bins):
        stop_bin2  = bin1+CORR_BINS+1 if bin1+CORR_BINS+1<N_bins else N_bins
        for bin2 in range(bin1,stop_bin2):
            # print bin1, bin2
            P_der_matrix[bin1,bin2]=DER(k,mu,var_num,bin1,bin2)
            P_der_matrix[bin2,bin1]=P_der_matrix[bin1,bin2]


# Compute Trace:
cdef double trace(double k, double mu, int var1, int var2):
    cdef double trace = 0
    cdef np.intp_t a, b, c, d
    if "correlations" in typeFM:
        # tick = time.clock()
        inverse_C = inverse_matrix_C(k, mu) * sqrt_volume_shells
        # inverse_C = inverse_matrix_C(k, mu)
        derivative_matrices(k, mu, var1, P_der_1_v)
        derivative_matrices(k, mu, var2, P_der_2_v)
        # tock = time.clock()
        # print "Inverse and derivatives: %g sec" %(tock-tick)

        # Numpy trace:
        # tick = time.clock()
        trace = np.dot(P_der_1, np.dot(inverse_C, np.dot(P_der_2,inverse_C))).trace()
        # tock = time.clock()
        # print "Numpy trace: %g sec --> %g" %(tock-tick, trace)

        # # Optimized Cython trace: (optimized con cavolo...)
        # inverse_C_v = inverse_matrix_C(k, mu)
        # tick = time.clock()
        # trace = 0
        # for a in range(N_bins):
        #     for b in range(N_bins):
        #         for c in range(N_bins):
        #             for d in range(N_bins):
        #                 trace=trace + P_der_1_v[a,b]*inverse_C_v[b,c]*P_der_2_v[c,d]*inverse_C_v[d,a] * sqrt( sqrt( vol_shell(a)*vol_shell(b)*vol_shell(c)*vol_shell(d)) )
        # tock = time.clock()
        # print "Cython trace: %g sec --> %g" %(tock-tick,trace)
        return trace
    #
    # Without correlations there is no Trace to compute,
    # just a sum over bins:
    #
    else:
        # Check if some element of the matrix is zero anyway:
        if var1>=N_cosm_vars_max and var2>=N_cosm_vars_max and var1!=var2:
            return 0.
        else:
            # Optimise the sum over bins:
            if var1>N_cosm_vars_max:
                bin = var1-N_cosm_vars_max # bin_bias
                # return DER(k, mu, var1, bin, bin) * DER(k, mu, var2, bin, bin) * (n_dens_c[bin]/(n_dens_c[bin]*observed_spectrum(bin,bin,k,mu)+1))**2
                return vol_shell(bin) * DER(k, mu, var1, bin, bin) * DER(k, mu, var2, bin, bin) * (n_dens_c[bin]/(n_dens_c[bin]*observed_spectrum(bin,bin,k,mu)+1))**2
            elif var2>N_cosm_vars_max:
                bin = var2-N_cosm_vars_max # bin_bias
                # return  DER(k, mu, var1, bin, bin) * DER(k, mu, var2, bin, bin) * (n_dens_c[bin]/(n_dens_c[bin]*observed_spectrum(bin,bin,k,mu)+1))**2
                return vol_shell(bin) * DER(k, mu, var1, bin, bin) * DER(k, mu, var2, bin, bin) * (n_dens_c[bin]/(n_dens_c[bin]*observed_spectrum(bin,bin,k,mu)+1))**2
            else:
                # tick = time.clock()
                result = 0.
                for bin in range(N_bins):
                    # result += DER(k, mu, var1, bin, bin) * DER(k, mu, var2, bin, bin) * (n_dens_c[bin]/(n_dens_c[bin]*observed_spectrum(bin,bin,k,mu)+1))**2
                    result += vol_shell(bin) * DER(k, mu, var1, bin, bin) * DER(k, mu, var2, bin, bin) * (n_dens_c[bin]/(n_dens_c[bin]*observed_spectrum(bin,bin,k,mu)+1))**2
                # tock = time.clock()
                # print "Everything: %g sec" %(tock-tick)
                return result

def trace_py(k,mu,var1,var2):
    return trace(k,mu,var1,var2)


# ------------------------------------------
# Adding redshift dependence:
# ------------------------------------------
cdef:
    double[::1] Growth_der_redshift # N_bins


def compute_data_adding_z():
    global Growth_der_redshift
    Growth_der_redshift = Fun_der_redshift(Growth,z_avg)


def obsSpectr_redshift_der(bin1,bin2,k,mu,der_bin):
    result = 0.
    if bin1==der_bin:
        result+=observed_spectrum(bin1,bin2,k,mu) * (Growth_bins[bin1] * Growth_der_redshift[bin1])
    if bin2==der_bin:
        result+=observed_spectrum(bin1,bin2,k,mu) * (Growth_bins[bin2] * Growth_der_redshift[bin2])
    return result

def sigma_redshift(bin):
    return (z_in[bin+1]-z_in[bin])/2.


def trace_adding_z(double k, double mu, int var1, int var2):
    N_tot_vars_max = N_cosm_vars_max + N_bins

    # Check if both are cosmo. params:
    if var1<N_tot_vars_max and var2<N_tot_vars_max:
        return trace(k,mu,var1,var2)

    inverse_C = inverse_matrix_C(k, mu) * sqrt_volume_shells
    global P_der_1, P_der_2, P_der_1_v, P_der_2_v

    # Compute first derivative matrix:
    if var1>=N_tot_vars_max:
        redshift1 = var1 - N_tot_vars_max
        P_der_1 = np.zeros((N_bins,N_bins))
        P_der_1_v = P_der_1
        for bin2 in range(N_bins):
            if bin2!=redshift1: # zero w/o correlations
                P_der_1_v[redshift1,bin2] = obsSpectr_redshift_der(redshift1,bin2,k,mu,redshift1)
                P_der_1_v[bin2,redshift1] = P_der_1_v[redshift1,bin2]
            else:
                P_der_1_v[bin2,bin2] = 2 * obsSpectr_redshift_der(bin2,bin2,k,mu,bin2)
    else:
        derivative_matrices(k, mu, var1, P_der_1_v)

    # Compute second derivative matrix:
    if var2>=N_tot_vars_max:
        redshift2 = var2 - N_tot_vars_max
        P_der_2 = np.zeros((N_bins,N_bins))
        P_der_2_v = P_der_2
        for bin2 in range(N_bins):
            if redshift2!=bin2:
                P_der_2_v[redshift2,bin2] = obsSpectr_redshift_der(redshift2,bin2,k,mu,redshift2)
                P_der_2_v[bin2,redshift2] = P_der_2_v[redshift2,bin2]
            else:
                P_der_2_v[bin2,bin2] = 2 * obsSpectr_redshift_der(bin2,bin2,k,mu,bin2)
    else:
        derivative_matrices(k, mu, var2, P_der_2_v)

    #Compute trace:
    result = np.dot(P_der_1, np.dot(inverse_C, np.dot(P_der_2,inverse_C))).trace()

    return result

# ------------------------------------------
# ------------------------------------------


# # Just for the adaptive sampling and interpolation: (ks is a vector!)
# # It computes the log for a better sampling
# def trace_log_args(ks,**args):
#     Nk = ks.shape[0]
#     args_names = ['mu', 'var1', 'var2']
#     args_values = [0.]*3
#     for i, arg_name in enumerate(args_names):
#         if arg_name in args:
#             args_values[i] = args[arg_name]
#         else:
#             print "Error trace: not all the args received!"
#     results = np.empty(Nk)
#     for i, k in enumerate(ks):
#         result = trace(k,args_values[0],args_values[1],args_values[2])
#         results[i] = result
#         # if result>0:
#         #     results[i] = np.log10(result)
#         # else:
#         #     results[i] = -np.log10(-result)
#     return results

# cdef double trace_part(double k, double mu, int var1, int var2, int bin_kmax):
#     inverse_C_v = inverse_matrix_C(k, mu)
#     derivative_matrices(k, mu, var1, P_der_1_v)
#     derivative_matrices(k, mu, var2, P_der_2_v)

#     # Optimized Cython trace:
#     cdef double trace = 0
#     cdef np.intp_t a, b, c, d
#     for a in range(0,bin_kmax+1):
#         for b in range(0,bin_kmax+1):
#             for c in range(0,bin_kmax+1):
#                 for d in range(0,bin_kmax+1):
#                     if ( a==bin_kmax or b==bin_kmax or c==bin_kmax or d==bin_kmax ):
#                         trace=trace + P_der_1_v[a,b]*inverse_C_v[b,c]*P_der_2_v[c,d]*inverse_C_v[d,a] * sqrt( sqrt( vol_shell_original(a)*vol_shell_original(b)*vol_shell_original(c)*vol_shell_original(d)) )
#     return(trace)


cdef double trace_part(double k, double mu, int var1, int var2, int bin_kmax):
    inverse_C_v = inverse_matrix_C(k, mu)
    derivative_matrices(k, mu, var1, P_der_1_v)
    derivative_matrices(k, mu, var2, P_der_2_v)

    # Optimized Cython trace:
    cdef double trace = 0
    cdef np.intp_t a, b, c, d
    for a in range(bin_kmax,N_bins):
        for b in range(bin_kmax,N_bins):
            for c in range(bin_kmax,N_bins):
                for d in range(bin_kmax,N_bins):
                    if ( a==bin_kmax or b==bin_kmax or c==bin_kmax or d==bin_kmax ):
                        trace=trace + P_der_1_v[a,b]*inverse_C_v[b,c]*P_der_2_v[c,d]*inverse_C_v[d,a] * sqrt( sqrt( vol_shell(a)*vol_shell(b)*vol_shell(c)*vol_shell(d)) )
    return(trace)


# #--------------------------------------------------------------
# #--------------------------------------------------------------
# # Trace interpolation: (just for performance opt.)
# #--------------------------------------------------------------
# #--------------------------------------------------------------
# cdef enum:
#     max_N_vars = 50
# cdef:
#     interpolation_tools_2D trace_tools[max_N_vars][max_N_vars]


# def trace_array(vect_k, vect_mu, var1, var2):
#     cdef:
#          int Nk = vect_k.shape[0], Nmu = vect_mu.shape[0]
#     results = np.empty((Nk,Nmu))
#     cdef:
#         double[::1] vect_k_c = vect_k
#         double[::1] vect_mu_c = vect_mu
#         double[:,::1] results_c = results
#     for ik in range(Nk):
#         for imu in range(Nmu):
#             results_c[ik,imu] = trace(vect_k_c[ik],vect_mu_c[imu],var1,var2)
#     return results


# def init_Trace(Nk=50,Nmu=4):
#     print "\nComputing, storing and interpolating traces:"

#     # Read content traces folder:
#     OUTPUT_PATH_TRACE = "INPUT/traces/"
#     for (dirpath, dirnames, filenames) in walk(OUTPUT_PATH_TRACE):
#         names = [ fi for fi in filenames if fi.endswith(".csv") ]
#         break

#     vect_k = np.linspace(1e-3,0.2,Nk)
#     vect_mu = np.linspace(-1.,1.,Nmu)
#     index = 0
#     for var1 in range(N_tot_vars):
#         for var2 in range(var1,N_tot_vars):
#             file_name = "Nk%d_Nmu%d_vars_%d-%d.csv" %(Nk,Nmu,var1,var2)
#             if file_name in names:
#                 if index==0:
#                     print "Traces computed previously. Importing from file..."
#                 flat_arr = np.loadtxt(open(OUTPUT_PATH_TRACE+file_name,"rb")).flatten()
#             else:
#                 start = time.time()
#                 flat_arr = trace_array(vect_k,vect_mu,var1,var2).T.flatten()
#                 np.savetxt(OUTPUT_PATH_TRACE+"Nk%d_Nmu%d_vars_%d-%d.csv" %(Nk,Nmu,var1,var2), flat_arr)
#                 run_time = time.time()-start
#                 remaining = run_time* (N_tot_vars*(N_tot_vars+1)/2 - index - 1) /60.
#                 print "- vars %d and %d: %g sec. \t --> ~%.0f min %d sec to go" %(var1,var2,run_time,remaining,int((remaining -int(remaining)) *60) )
#             alloc_interp_GSL_2D(vect_k, vect_mu, flat_arr, &trace_tools[var1][var2])
#             index+=1
#     print "Done!"

# def prova(x,y):
#     x_max, y_max = x[-1], y[0]
#     y = y/y_max * x_max
#     return y*y*y*y

# def adapt_trace_term(var1, var2, Nk_start=50, mu=-1., tol=0.05, min_points=10, max_level=20):
#     args = {'mu':mu, 'var1':var1, 'var2':var2}
#     vect_k = np.linspace(1e-3,0.2,Nk_start)
#     return adapt_sampl.sample_function(trace_log_args, vect_k, tol, min_points, max_level, prova, **args)

# def adapt_trace_tot(var1, var2, Nk_start=15, Nmu=10, tol=0.05, min_points=10, max_level=20):
#     # Compute optimised k for mu = -1.0:
#     samples_k, logTr_mu1 = adapt_trace_term(var1, var2, Nk_start, -1., tol, min_points, max_level)
#     print "Numer samples: %d" %(samples_k.shape[0])

#     # Compute trace samples for other values of mu:
#     Nk_sample = samples_k.shape[0]
#     logTr = np.empty((Nk_sample,Nmu))
#     logTr[:,0] = logTr_mu1.T
#     vect_mu = np.linspace(-1.,1.,Nmu)
#     for i, mu in enumerate(vect_mu[1:]):
#         args = {'mu':mu, 'var1':var1, 'var2':var2}
#         logTr[:,i+1] = trace_log_args(samples_k,**args).T
#         print ".",

#     # Interpolate k-direction with Akima1DInterpolator: (avoid spline mess)
#     print "\nAkima... ",
#     np.savetxt("temp.csv",logTr)
#     Nk = 3000
#     vect_k = np.linspace(1e-3,0.2,Nk)
#     final_trSamples = np.empty((Nk,Nmu))
#     for i, mu in enumerate(vect_mu):
#         final_trSamples[:,i] = Akima1DInterpolator(samples_k, logTr[:,i].T)(vect_k)
#     # Revert log:
#     # idxs1, idxs2 = (final_trSamples >= 0.).nonzero(), (final_trSamples < 0.).nonzero()
#     # final_trSamples[idxs1] = np.power(10,final_trSamples[idxs1])
#     # final_trSamples[idxs2] = - np.power(10,-final_trSamples[idxs2])

#     # Finally interpolate everything with GSL 2D: (spline cubic)
#     print "GSL..."
#     alloc_interp_GSL_2D(vect_k, vect_mu, final_trSamples.T.flatten(), &trace_tools[var1][var2])


# def init_Trace_term(var1, var2, Nk=50, Nmu=4):
#     # Read content traces folder:
#     OUTPUT_PATH_TRACE = "INPUT/traces/"
#     for (dirpath, dirnames, filenames) in walk(OUTPUT_PATH_TRACE):
#         names = [ fi for fi in filenames if fi.endswith(".csv") ]
#         break

#     vect_k = np.linspace(1e-3,0.2,Nk)
#     vect_mu = np.linspace(-1.,1.,Nmu)

#     file_name = "Nk%d_Nmu%d_vars_%d-%d.csv" %(Nk,Nmu,var1,var2)
#     if file_name in names:
#         print "Trace interpolated previously. Importing from file..."
#         flat_arr = np.loadtxt(open(OUTPUT_PATH_TRACE+file_name,"rb")).flatten()
#     else:
#         start = time.time()
#         flat_arr = trace_array(vect_k,vect_mu,var1,var2).T.flatten()
#         np.savetxt(OUTPUT_PATH_TRACE+"Nk%d_Nmu%d_vars_%d-%d.csv" %(Nk,Nmu,var1,var2), flat_arr)
#         run_time = time.time()-start
#         print "Computing trace: %g sec." %(run_time)
#     alloc_interp_GSL_2D(vect_k, vect_mu, flat_arr, &trace_tools[var1][var2])


# cdef double interp_trace(double k, double mu, int var1, int var2):
#     return eval_interp_GSL_2D(k, mu, &trace_tools[var1][var2])

# def interp_trace_py(k,mu,var1,var2):
#     return interp_trace(k,mu,var1,var2)

#--------------------------------------------------------------
#--------------------------------------------------------------
# Computation of the Fisher Matrix: (with GSL)
#--------------------------------------------------------------
#--------------------------------------------------------------

# Default values are set in set_kargs_FM():
rel_prec, abs_prec, rel_prec_redshift, abs_prec_redshift, rel_prec_redshift_cross, abs_prec_redshift_cross = 0., 0., 0., 0., 0., 0.


# Integration variables:
cdef enum:
    max_alloc_const = 500
    max_alloc_const_K = 10000
cdef:
    size_t MAX_ALLOC = max_alloc_const
    size_t MAX_ALLOC_K = max_alloc_const_K


cdef:
    gsl_integration_workspace * W_k
    gsl_integration_workspace * W_mu
W_k = gsl_integration_workspace_alloc(MAX_ALLOC_K)
W_mu = gsl_integration_workspace_alloc(MAX_ALLOC)

# Integral arguments:
cdef double argument_mu(double mu, void *input): #k, var1, var2, #bin_kmax
    cdef:
        double *params = <double*> input
        double k = params[0]
    # OBSOLETE:
    if abs(mode_kmax)<1e-10: # is zero
        return( k*k * trace_part(k,mu,<int>params[1],<int>params[2],<int>params[3]) )
    else:
        # if interpolate_Trace:
        #     return k*k * interp_trace(k,mu,<int>params[1],<int>params[2])
        # else:
        return k*k * trace_adding_z(k,mu,<int>params[1],<int>params[2])

cdef double argument_k(double k, void *input): #var1, var2, bin_kmax
    cdef double params[4]
    cdef double *vars = <double*>input
    var1, var2 = vars[0], vars[1]
    cdef gsl_function F_mu
    F_mu.function = &argument_mu
    params[0], params[1], params[2], params[3] = k, var1, var2, vars[2]

    # Redshift dependence: (adjust integral prec.)
    FM_mu_absPrec, FM_mu_relPrec = abs_prec, rel_prec
    N_tot_vars_max = N_cosm_vars_max + N_bins
    if var1>=N_tot_vars_max and var2>=N_tot_vars_max:
        FM_mu_absPrec, FM_mu_relPrec = abs_prec_redshift, rel_prec_redshift
    elif var1>=N_tot_vars_max or var2>=N_tot_vars_max:
        FM_k_absPrec, FM_k_relPrec = abs_prec_redshift_cross, rel_prec_redshift_cross

    cdef double result = eval_integration_GSL(-1., 1., FM_mu_absPrec, FM_mu_relPrec, params, W_mu, &F_mu, MAX_ALLOC)
    return result




#------------------------
# FISHER MATRIX element:
#------------------------
cdef:
    double mode_kmax # 0 for k_max(z), number for fixed k_max

AP_flag = False
zDist_flag = False
interpolate_Trace = False
# typeFM = "correlations+windowFun"
typeFM = "uncorr"

k_min_hard = 5e-3
k_max_hard = 0.5


def set_kargs_FM(**kargs):
    global AP_flag, typeFM, CORR_BINS, mode_kmax, zDist_flag

    # AP effect:
    AP_flag = kargs["check_AP"] if "check_AP" in kargs else False
    zDist_flag = kargs["zDistortion"] if "zDistortion" in kargs else False


    # Type of FM:
    typeFM = kargs["typeFM"] if "typeFM" in kargs else "uncorr"

    # Adjustable k_max FM: (not-working at the moment..)
    mode_kmax = kargs["fixed_kmax"] if "fixed_kmax" in kargs else 0.2

    # Number of correlated bins:
    if "num_corr_bins" in kargs:
        CORR_BINS = kargs["num_corr_bins"]
    elif "Delta_z_correlations" in kargs:
        delta_bin = z_in[1]-z_in[0]
        division = kargs["Delta_z_correlations"] / delta_bin
        if division>1.:
            CORR_BINS = int(division)+1 if (division%int(division))>0.25 else int(division)
        else:
            CORR_BINS = 1 # Keep at least one correlation
    else:
        CORR_BINS = N_bins
    kargs["num_corr_bins"] = CORR_BINS

    # k_min and k_max:
    global k_min_hard, k_max_hard
    k_min_hard = kargs["k_min_hard"] if "k_min_hard" in kargs else 5e-3
    k_max_hard = kargs["k_max_hard"] if "k_max_hard" in kargs else 0.5

    # Integral's precisions: (add separate k and mu...)
    global rel_prec, abs_prec
    global rel_prec_redshift, abs_prec_redshift, rel_prec_redshift_cross, abs_prec_redshift_cross
    rel_prec = kargs["rel_prec"] if "rel_prec" in kargs else 1e-2
    abs_prec = kargs["abs_prec"] if "abs_prec" in kargs else 100.
    rel_prec_redshift = kargs["rel_prec_redshift"] if "rel_prec_redshift" in kargs else 1e-4
    abs_prec_redshift = kargs["abs_prec_redshift"] if "abs_prec_redshift" in kargs else 1e-2
    rel_prec_redshift_cross = kargs["rel_prec_redshift_cross"] if "rel_prec_redshift_cross" in kargs else 1e-2
    abs_prec_redshift_cross = kargs["abs_prec_redshift_cross"] if "abs_prec_redshift_cross" in kargs else 1



    # Redshift dependence:
    if "z_dependence" not in kargs:
        kargs["z_dependence"] = False
    elif kargs["z_dependence"]==True:
        # Compute only at the beginning, not for every FM element:
        if "FM_elem" not in kargs:
            compute_data_adding_z()

    return kargs



def fisher_matrix_element(int var1, int var2, **kargs):
    kargs["FM_elem"] = True
    kargs = set_kargs_FM(**kargs)

    cdef gsl_function F_k
    F_k.function = &argument_k
    cdef:
        double params[3]
        double FM_elem = 0
    params[0], params[1] = var1, var2
    FM_k_absPrec, FM_k_relPrec = abs_prec, rel_prec

    # Adding redshift dependence:
    N_tot_vars_max = N_cosm_vars_max + N_bins
    if var1>=N_tot_vars_max and var2>=N_tot_vars_max:
        FM_k_absPrec, FM_k_relPrec = abs_prec_redshift, rel_prec_redshift
        if var1==var2:
            FM_elem += 1./sigma_redshift(var1-N_tot_vars_max)**2
    elif var1>=N_tot_vars_max or var2>=N_tot_vars_max:
        FM_k_absPrec, FM_k_relPrec = abs_prec_redshift_cross, rel_prec_redshift_cross


    # Sum all the integrals with different k_max: (NOT WORKING)
    # trace_part() takes care of selecting the appropriate terms of the trace
    if abs(mode_kmax)<1e-10: # is zero
        for bin_max in range(N_bins):
            if (k_max_data[bin_max]<k_max_hard):
                k_max_int=k_max_data[bin_max]
            else:
                k_max_int=k_max_hard
            params[2] = bin_max
            FM_elem += 1./(8*np.pi**2) * eval_integration_GSL(k_min_hard, k_max_int, FM_k_absPrec, FM_k_relPrec, params, W_k, &F_k, MAX_ALLOC_K)
    else:
        FM_elem += 1./(8*np.pi**2) * eval_integration_GSL(k_min_hard, mode_kmax, FM_k_absPrec, FM_k_relPrec, params, W_k, &F_k, MAX_ALLOC_K)

    return FM_elem



#------------------------
# Computation of FM:
#------------------------
#
# Types available:
#  - uncorrelated
#  - windowFun
#  - correlations+windowFun
#  - correlations (not really...)
#


def FM(FMname=time.strftime("%d-%m-%y")+"_"+time.strftime("%H-%M"),**kargs):

    kargs = set_kargs_FM(**kargs)
    print "\nComputing Fisher matrix... (type: %s)" %(typeFM)
    if kargs["z_dependence"]==True:
        print "--> Redshift dependence added to FM"
    if "correlations" in typeFM:
        print "--> Num. correlated bins: %d" %(CORR_BINS)



    FM_dim_max = N_cosm_vars_max + N_bins if not kargs["z_dependence"] else N_cosm_vars_max + 2*N_bins
    FM_dim = N_cosm_vars+N_bins if not kargs["z_dependence"] else N_cosm_vars+2*N_bins
    FM = np.zeros((FM_dim,FM_dim))

    # Resetting matrices:
    global P_der_1, P_der_1_v, P_der_2, P_der_2_v, C, C_v
    P_der_1, P_der_2 = np.zeros([N_bins, N_bins]), np.zeros([N_bins, N_bins])
    P_der_1_v, P_der_2_v = P_der_1, P_der_2
    C = np.zeros([N_bins, N_bins])
    C_v = C

    # Checking AP and import AP_integrals:
    if AP_flag:
        if "correlations" in typeFM and (CORR_BINS>1 or N_bins!=14):
            print "Error: correlations not computed for the selected case. Use the AP module before."
            return -1
        print "--> AP effect included"
        import_AP_int()

    # Some log:
    if "log_file" in kargs:
        logFile = open(kargs["log_file"],'a')
        logFile.write("\nComputing Fisher matrix... (type: %s)\n" %(typeFM))
        if kargs["z_dependence"]==True:
            logFile.write("--> Redshift dependence added to FM\n")
        if "correlations" in typeFM:
            logFile.write("--> Num. correlated bins: %d\n" %(CORR_BINS) )
        if AP_flag:
            logFile.write("--> AP effect included\n")


    # Computing:
    var_numbers = FM_vars_numbers + range(N_cosm_vars_max,FM_dim_max)
    for i, var1 in enumerate(var_numbers):
        for j, var2 in zip(range(i,FM_dim),var_numbers[i:]):
            tick = time.time()
            FM[i,j]=fisher_matrix_element(var1,var2,**kargs)
            tock = time.time()
            FM[j,i]=FM[i,j]
            if "log_file" in kargs:
                logFile.write("\n(%d,%d) --> %g sec.\n" %(i,j,tock-tick))
            if "verbose" in kargs:
                print "(%d,%d) --> %g sec." %(i,j,tock-tick)
            if 'correlations' in typeFM:
                np.savetxt("OUTPUT/FM_%s_AP%d-zDep%d-%dbins-%s.csv" %(FMname,AP_flag,kargs["z_dependence"],N_bins,typeFM), FM)
                if N_bins>25:
                    print "(%d,%d) --> %g sec." %(i,j,tock-tick)
    if 'correlations' not in typeFM:
        np.savetxt("OUTPUT/FM_%s_AP%d-zDep%d-%dbins-%s.csv" %(FMname,AP_flag,kargs["z_dependence"],N_bins,typeFM), FM)
    return FM



def set_typeFM(type_FM_input):
    global typeFM
    typeFM = type_FM_input


# def FM_plusRedshift(FMname="test",**kargs):

#     set_kargs_FM(**kargs)

#     print "\nComputing Fisher matrix (redshift). Type: %s" %(typeFM)
#     N_tot_vars = N_cosm_vars + N_bins
#     FM_dimension = N_tot_vars + N_bins
#     FM = np.zeros([FM_dimension,FM_dimension])

#     # Resetting matrices:
#     global P_der_1, P_der_1_v, P_der_2, P_der_2_v, C, C_v
#     P_der_1, P_der_2 = np.zeros([N_bins, N_bins]), np.zeros([N_bins, N_bins])
#     P_der_1_v, P_der_2_v = P_der_1, P_der_2
#     C = np.zeros([N_bins, N_bins])
#     C_v = C

#     if Delta_z_correlations<0:
#         num_corr_bins = N_bins
#     else:
#         delta_bin = z_in[1]-z_in[0]
#         division = Delta_z_correlations / delta_bin
#         if division>1.:
#             num_corr_bins = int(division)+1 if (division%int(division))>0.25 else int(division)
#         else:
#             num_corr_bins = 1
#         print "Num. correlated bins: %d" %(num_corr_bins)


#     # Computing:
#     var_numbers = FM_vars_numbers + range(N_cosm_vars_max,N_cosm_vars_max+2*N_bins)
#     for i, var1 in enumerate(var_numbers):
#         for j, var2 in zip(range(i,FM_dimension),var_numbers[i:]):
#             FM[i,j]=fisher_matrix_element(var1,var2,check_AP,type_FM_input,num_corr_bins,fixed_kmax)
#             FM[j,i]=FM[i,j]
#             if 'correlations' in type_FM_input:
#                 np.savetxt("OUTPUT/FM_plusRedshift_%s_AP%d-%dbins-%s.csv" %(FMname,check_AP,N_bins,type_FM_input), FM)
#     if 'correlations' not in type_FM_input:
#         np.savetxt("OUTPUT/FM_plusRedshift_%s_AP%d-%dbins-%s.csv" %(FMname,check_AP,N_bins,type_FM_input), FM)
#     return FM





# include "libraries/adding-z.pyx"


# #################
# PLOTTING TRACE:
# #################

# def plot_trace(int var1, int var2, N_k=100,N_mu=10,k_min=1e-3,k_max=0.5):
#     k_vect = np.linspace(k_min,k_max,N_k)
#     # k_vect = np.logspace(np.log10(k_min),np.log10(k_max),N_k)
#     mu_vect = np.linspace(-1.,1.,N_mu)
#     samples = np.zeros([N_k,N_mu])
#     cdef:
#         double[:,::1] samples_c = samples
#         double[::1] k_vect_c = k_vect
#         double[::1] mu_vect_c = mu_vect
#     time_count = 0
#     total_start = time.time()
#     for i_k in range(N_k):
#         for i_mu in range(N_mu):
#             samples_c[i_k,i_mu] = trace(k_vect_c[i_k],mu_vect_c[i_mu],var1,var2)
#     total_stop = time.time()
#     print "Trace computation: %g seconds" %(total_stop-total_start)
#     np.savetxt("OUTPUT/trace/trace_3D_vars_%d%d.csv" %(var1,var2),samples)

#     # ######################
#     # # Plot this thing...
#     # # Result ---> it's a mess ;D
#     # ######################
#     # from mpl_toolkits.mplot3d import Axes3D
#     # from matplotlib import cm
#     # from matplotlib.ticker import LinearLocator, FormatStrFormatter
#     # fig = pl.figure()
#     # ax = fig.gca(projection='3d')
#     # # X = np.arange(-5, 5, 0.25)
#     # # Y = np.arange(-5, 5, 0.25)
#     # # R = np.sqrt(X**2 + Y**2)
#     # # Z = np.sin(R)
#     # mu_vect_m, k_vect_m  = np.meshgrid(mu_vect, k_vect)
#     # surf = ax.plot_surface(k_vect_m, mu_vect_m, np.log10(samples), rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
#     # # # ax.set_zlim(-1.01, 1.01)
#     # ax.zaxis.set_major_locator(LinearLocator(10))
#     # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#     # # ax.zaxis.set_scale('log')
#     # fig.colorbar(surf, shrink=0.5, aspect=5)
#     # fig.savefig('plots/trace/trace_3D_vars_%d%d.pdf'%(var1,var2))

#     # fig2=pl.figure()
#     # ax1=fig2.add_subplot(111)
#     # ax1.plot(k_vect,samples[:,0],'r-',label="Trace along $\\mu=-1$")
#     # # ax1.plot(vect_k,class_fct['P_0'](vect_k),'b-',label="Class Spectrum")
#     # ax1.grid(True)
#     # ax1.legend(loc='best')
#     # ax1.set_yscale('log')
#     # # ax1.set_xlabel("$k$ [$h$/Mpc]")
#     # # ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
#     # fig2.savefig('plots/trace/trace_along_mu_vars%d%d.pdf'%(var1,var2))
#     # np.savetxt("OUTPUT/trace/trace_along_mu_vars%d%d.csv"%(var1,var2),np.column_stack((k_vect,samples[:,0])))
#     return k_vect, samples[:,0]



# include "libraries/plot_stuff.pxi"

# # ###########################
# # Plot k-derivative comparison:
# # ###########################
# def plot_k_der_compare(int var_num, bin1, bin2, double k_min=0.001, double k_max=0.5, int N_k=1000):
#     # N_k, N_mu = 1000, 10
#     k_vect = np.linspace(k_min,k_max,N_k)
#     samples = np.empty(N_k)
#     samples2 = np.empty(N_k)
#     # F_k.function = &argument_k
#     # F_mu.function = &argument_mu
#     # sdfsdf[0] = 3.
#     mu = 1.
#     for i_k in range(N_k):
#         samples[i_k] = windowed_DER_k(k_vect[i_k],mu,bin1,bin2,var_num)
#         samples2[i_k] = zero_spectrum_der_k(k_vect[i_k])*k_der(mu, k_vect[i_k], bin1, var_num)
#     fig2=pl.figure()
#     ax1=fig2.add_subplot(111)
#     ax1.plot(k_vect,samples,'r-',label="Windowed")
#     ax1.plot(k_vect,samples2,'b-',label="Normal")
#     ax1.grid(True)
#     ax1.legend(loc='best')
#     # ax1.set_yscale('log')
#     # ax1.set_xlabel("$k$ [$h$/Mpc]")
#     # ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
#     fig2.savefig('plots/compare_derivatives_bins%d%d.pdf' %(bin1,bin2))
#     # np.savetxt("output/FM_arg_k_integral_vars%d%d.csv" %(var1,var2),np.column_stack((k_vect,samples)))
#     return

