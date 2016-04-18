###########
# Brunch:
###########
# - modified Trace that set k_max for the lower index (commented there is also the higher case)



# Importing modules and defining default data:

include "libraries/header.pxi"

include "libraries/import_CLASS.pxi"

include "libraries/analytical_fcts.pxi"


#**************************************************************
# Importing windowed data and interpolation: (with GSL)
#**************************************************************

# Tools for GLS interpolation:
cdef enum:
    max_N_bins = 50

cdef:
    interpolation_tools integral_1[max_N_bins][max_N_bins][N_vars]
    interpolation_tools integral_DER[max_N_bins][max_N_bins]

# Read file names:
PATH_INPUT_WIND_SPECTRA = "INPUT/modified_spectra/"
def import_wind_spectra():
    # global k_min, k_max
    for (dirpath, dirnames, filenames) in walk(PATH_INPUT_WIND_SPECTRA):
        names = [ fi for fi in filenames if fi.endswith(".csv") ]
        break
    # Importing files:
    print "Importing and interpolating windowed spectra. . .",
    j, count1, countDER = 0, 0, 0
    for filename in names:
        # Some log:
        if j%10==0:
            print".",
        j+=1
        file_data = np.loadtxt(open(PATH_INPUT_WIND_SPECTRA+filename,"rb"),delimiter=" ")
        k_data = np.array(file_data[:,0])
        # Determine which kind of data there is in the file:
        if filename.endswith("int1.csv"): # (var_num from 0 to 4)
            m = match(r"([0-9]*)_([0-9]*)-var([0-4])_int1", filename)
            bin1, bin2, var_num = int(m.group(1)), int(m.group(2)), int(m.group(3))
            data_y = np.array(file_data[:,1])
            alloc_interp_GSL( k_data, data_y, &integral_1[bin1][bin2][var_num])
            count1+=1
        elif filename.endswith("intDER.csv"):
            m = match(r"([0-9]*)_([0-9]*)_intDER", filename)
            bin1, bin2 = int(m.group(1)), int(m.group(2))
            data_DER = np.array(file_data[:,1])
            # Complete missing data to k=k_max:
            added_k = np.linspace(k_data[-1],k_max,1001)[1:]
            added_zeros = np.zeros(1000)
            alloc_interp_GSL( np.concatenate((k_data,added_k)), np.concatenate((data_DER,added_zeros)), &integral_DER[bin1][bin2])
            countDER+=1
        # Save max_k and min_k:
        # if j==1:
        #     k_min, k_max = k_data[0], k_data[-1]
    N1_expected = N_bins*(N_bins+1)/2*(5)
    NDer_expected = N_bins*N_bins
    if ( count1!=N1_expected or countDER!=NDer_expected ):
        print "\n\n !!!!  ---> Some windowed spectra missing: int1=%d, int23=%d --> Expected: (%d, %d)  <---  !!!!" %(count1,countDER,N1_expected,NDer_expected)
    else:
        print "Successfull!"



# Initialise everything:
def all():
    import_CLASS_data()
    compute_survey_DATA()
    import_wind_spectra()


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
cdef double windowed_Spectrum(double k, int bin1, int bin2):
    if bin1<=bin2:
        return eval_interp_GSL(k, &integral_1[bin1][bin2][0])
    else:
        return eval_interp_GSL(k, &integral_1[bin2][bin1][0])

cdef double windowed_DER(double k, int bin1, int bin2, int var): #var [0-3]
    if bin1<=bin2:
        return eval_interp_GSL(k, &integral_1[bin1][bin2][var+1])
    else:
        return eval_interp_GSL(k, &integral_1[bin2][bin1][var+1])

cdef double windowed_DER_k(double k, double mu, int bin1, int bin2, int var): #var [3-5]
    W2_derW1 = eval_interp_GSL(k, &integral_DER[bin1][bin2])
    W1_derW2 = eval_interp_GSL(k, &integral_DER[bin2][bin1])
    return sqrt(vol_shell_original(bin1)*vol_shell_original(bin2)) / (4*PI**2) * ( k*mu*mu*(W2_derW1*(lnH_der_data[var][bin1]+lnD_der_data[var][bin1]) +W1_derW2*(lnH_der_data[var][bin2]+lnD_der_data[var][bin2])) -k*(W2_derW1*lnD_der_data[var][bin1] +W1_derW2*lnD_der_data[var][bin2] ) )


#--------------------------------------------------------------
# Contructing the 16 final derivatives for the Fisher Matrix:
#--------------------------------------------------------------
# Observed spectrum: (optimized!)
cdef double observed_spectrum(int bin1, int bin2, double k, double mu):
    return( Growth_bins[bin1]*Growth_bins[bin2] * bias_bins[bin1]*bias_bins[bin2]* (1+beta_bins[bin1]*mu**2)*(1+beta_bins[bin2]*mu**2) * windowed_Spectrum(k,bin1,bin2) )

# Observed terms: (optimized!)
# to avoid division by zero given by windowed_Spectrum with i!=j
cdef double observed_terms(int bin1, int bin2, double k, double mu):
    return( Growth_bins[bin1]*Growth_bins[bin2] * bias_bins[bin1]*bias_bins[bin2]* (1+beta_bins[bin1]*mu**2)*(1+beta_bins[bin2]*mu**2) )

# h, Omega_b and n_s: (optimized!) (var 0-2)
cdef double der_type_A(int bin1, int bin2, double k,double mu, int var_num):
    return  observed_terms(bin1, bin2, k, mu)*windowed_DER(k, bin1, bin2, var_num)

# Omega_m, w_p and w_1: (optimized!) (var 3-5)
cdef double der_type_B(int bin1, int bin2, double k, double mu, int var_num):
    cdef double CLASS_term
    if var_num==3: #Omega_m
        CLASS_term = windowed_DER(k, bin1, bin2, var_num)
    else:
        CLASS_term = 0
    cdef double beta_term = 1./(1+beta_bins[bin1]*mu**2)*( 2*beta_bins[bin1]*mu*mu_der(mu,bin1,var_num) + mu**2*Beta_der_data[var_num][bin1] ) +  1./(1+beta_bins[bin2]*mu**2)*( 2*beta_bins[bin2]*mu*mu_der(mu,bin2,var_num)+ mu**2*Beta_der_data[var_num][bin2] )
    cdef double derivative_k = windowed_DER_k(k,mu,bin1,bin2,var_num)


    # TEMPORARY REPLACEMENT FOR VOID and WRONG DATA in intDER:
    if (bin1==bin2 and k>0.43): #some instability after this point
        derivative_k =  zero_spectrum_der_k(k)*k_der(mu, k, bin1, var_num)


    # Pay attention to lnH_der_data that are computed in z_avg....!!!
    cdef np.intp_t avg_bin = (bin1+bin2)/2
    cdef double temp = observed_terms(bin1, bin2, k, mu) * ( CLASS_term + derivative_k*CHECK_DER_K) + observed_spectrum(bin1, bin2, k, mu) * ( lnG_der_data[var_num][bin1]+lnG_der_data[var_num][bin2] + lnH_der_data[var_num][avg_bin] - 2*lnD_der_data[var_num][avg_bin] + beta_term  )
    return(temp)

# Gamma: (optimized!) (var=6)
cdef double der_gamma(int bin1, int bin2, double k, double mu):
    # Pay attention to lnH_der_data that are computed in z_avg....
    return(observed_spectrum(bin1, bin2, k, mu) * (lnG_der_data[6][bin1]+lnG_der_data[6][bin2] + 1./(1+beta_bins[bin1]*mu**2)*(mu**2*Beta_der_data[6][bin1]) + 1./(1+beta_bins[bin2]*mu**2)*(mu**2*Beta_der_data[6][bin2])) )

# Bias: (9 derivatives) (bad optimized....) (var>=7)
cdef double der_bias(int bin1, int bin2, double k, double mu, int bin_bias) except -1:
    bias_term = lambda bin: 1/bias_bins[bin] - 1./(1+beta_bins[bin]*mu**2)*mu**2 * fnEv(Om_m_z_py,z=z_avg[bin],w_1=ref_values['w_1'],w_p=ref_values['w_p'],Om_m=ref_values['Om_m'])**ref_values['gamma'] /(bias_bins[bin]**2)
    if bin1==bin2 and bin1==bin_bias:
        return(observed_spectrum(bin1, bin2, k, mu) * 2*bias_term(bin_bias))
    elif bin1==bin_bias or bin2==bin_bias:
        return(observed_spectrum(bin1, bin2, k, mu) * bias_term(bin_bias))
    else:
        return(0.)


#--------------------------------------------------------------
# Constructing matrices and find the final Trace: (numpy+cython)
# (optimized with memory-views)
#--------------------------------------------------------------

# # Defining variables:
# P_der_1, P_der_2 = np.zeros([N_bins, N_bins]), np.zeros([N_bins, N_bins])
# C = np.zeros([N_bins, N_bins])

# # Noise matrix N:
# efficiency = 0.5 # Boh...
# n_dens = efficiency*np.array([3.56, 2.42, 1.81, 1.44, 0.99, 0.55, 0.29, 0.15])*1e-3
# inv_n_dens = 1./n_dens
# N = np.identity(N_bins) * inv_n_dens

# Compute inverse covariance matrix:
cdef double[:,::1] inverse_matrix_C(double k, double mu):
    cdef np.intp_t bin1, bin2
    for bin1 in range(N_bins):
        for bin2 in range(bin1,N_bins):
            C_v[bin1,bin2]= observed_spectrum(bin1, bin2, k, mu) + N_v[bin1,bin2]
            C_v[bin2,bin1]=C_v[bin1,bin2]
    return(np.linalg.inv(C))

# Compute matrices of derivatives: (var_num from 0 to 14)
cdef void derivative_matrices(double k, double mu, int var_num, double[:,::1] P_der_matrix):
    cdef np.intp_t bin1, bin2
    if var_num+1<=3: #h, Omega_b and n_s
        for bin1 in range(N_bins):
            for bin2 in range(bin1,N_bins):
                P_der_matrix[bin1,bin2]=der_type_A(bin1,bin2,k,mu,var_num)
                P_der_matrix[bin2,bin1]=P_der_matrix[bin1,bin2]
    elif var_num+1<=6: # Omega_m, w_p and w_1:
        for bin1 in range(N_bins):
            for bin2 in range(bin1,N_bins):
                P_der_matrix[bin1,bin2]=der_type_B(bin1,bin2,k,mu,var_num)
                P_der_matrix[bin2,bin1]=P_der_matrix[bin1,bin2]
    elif var_num+1==7: #gamma
        for bin1 in range(N_bins):
            for bin2 in range(bin1,N_bins):
                P_der_matrix[bin1,bin2]=der_gamma(bin1,bin2,k,mu)
                P_der_matrix[bin2,bin1]=P_der_matrix[bin1,bin2]
    else: #bias
        bin_bias = var_num-7
        for bin1 in range(N_bins):
            for bin2 in range(bin1,N_bins):
                P_der_matrix[bin1,bin2]=der_bias(bin1,bin2,k,mu,bin_bias)
                P_der_matrix[bin2,bin1]=P_der_matrix[bin1,bin2]
    return

# Compute Trace:
cdef double trace(double k, double mu, int var1, int var2):
    inverse_C_v = inverse_matrix_C(k, mu)
    derivative_matrices(k, mu, var1, P_der_1_v)
    derivative_matrices(k, mu, var2, P_der_2_v)

    # Optimized Cython trace:
    cdef double trace = 0
    cdef np.intp_t a, b, c, d
    for a in range(N_bins):
        for b in range(N_bins):
            for c in range(N_bins):
                for d in range(N_bins):
                    trace=trace + P_der_1_v[a,b]*inverse_C_v[b,c]*P_der_2_v[c,d]*inverse_C_v[d,a] * sqrt( sqrt( vol_shell_original(a)*vol_shell_original(b)*vol_shell_original(c)*vol_shell_original(d)) )
    return(trace)

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
                        trace=trace + P_der_1_v[a,b]*inverse_C_v[b,c]*P_der_2_v[c,d]*inverse_C_v[d,a] * sqrt( sqrt( vol_shell_original(a)*vol_shell_original(b)*vol_shell_original(c)*vol_shell_original(d)) )
    return(trace)


#--------------------------------------------------------------
#--------------------------------------------------------------
# Computation of the Fisher Matrix: (with GSL)
#--------------------------------------------------------------
#--------------------------------------------------------------

# Integration variables:
cdef enum:
    max_alloc_const = 500
    max_alloc_const_K = 10000
cdef:
    size_t MAX_ALLOC = max_alloc_const
    size_t MAX_ALLOC_K = max_alloc_const_K
    double rel_prec = 1e-2
    double abs_prec = 1e-6
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
    if abs(mode_kmax)<1e-10: # is zero
        return( k*k * trace_part(k,mu,<int>params[1],<int>params[2],<int>params[3]) )
    else:
        return( k*k * trace(k,mu,<int>params[1],<int>params[2]) )

cdef double argument_k(double k, void *input): #var1, var2, bin_kmax
    cdef double params[4]
    cdef double *vars = <double*>input
    cdef gsl_function F_mu
    F_mu.function = &argument_mu
    params[0], params[1], params[2], params[3] = k, vars[0], vars[1], vars[2]
    cdef double result = eval_integration_GSL(0., 1., abs_prec, rel_prec, params, W_mu, &F_mu, MAX_ALLOC)
    # print "(%g,   %g) "%(k, result)
    return 2*result

#------------------------
# FISHER MATRIX element:
#------------------------
cdef:
    double CHECK_DER_K = 0
    double mode_kmax # 0 for k_max(z), number for fixed k_max
k_min_hard = 0.001
k_max_hard = 0.5
def fisher_matrix_element(int var1, int var2, int check_K=0, double fixed_kmax=0):
    global CHECK_DER_K, mode_kmax
    CHECK_DER_K=check_K
    mode_kmax = fixed_kmax
    cdef gsl_function F_k
    F_k.function = &argument_k
    # Let's laugh... :/ --> IT WORKS!! :D (Apparently)
    cdef:
        double params[3]
        double FM_elem = 0
    params[0], params[1] = var1, var2
    # Sum all the integrals with different k_max:
    # trace_part() takes care of selecting the appropriate terms of the trace
    if abs(mode_kmax)<1e-10: # is zero
        for bin_max in range(N_bins):
            if (k_max_data[bin_max]<k_max_hard):
                k_max_int=k_max_data[bin_max]
            else:
                k_max_int=k_max_hard
            params[2] = bin_max
            FM_elem += 1./(8*np.pi**2) * eval_integration_GSL(k_min_hard, k_max_int, abs_prec, rel_prec, params, W_k, &F_k, MAX_ALLOC_K)
    else:
        FM_elem += 1./(8*np.pi**2) * eval_integration_GSL(k_min_hard, fixed_kmax, abs_prec, rel_prec, params, W_k, &F_k, MAX_ALLOC_K)
    return(FM_elem)


#------------------------
# Computation of FM:
#------------------------
def FM(int check_K=0, double fixed_kmax=0):
    FM = np.zeros([N_tot_vars,N_tot_vars])
    for var1 in range(N_tot_vars):
        for var2 in range(var1,N_tot_vars):
            start = time.clock()
            FM[var1,var2]=fisher_matrix_element(var1,var2,check_K,fixed_kmax)
            stop = time.clock()
            FM[var2,var1]=FM[var1,var2]
            np.savetxt("OUTPUT/FM_corr_K%d_kmax%g.csv" %(CHECK_DER_K,fixed_kmax), FM)
            print "(%d, %d) --> %g (%g sec.)" %(var1,var2,FM[var1,var2],(stop-start))
    return FM

def FM_var(int var1, int check_K=0):
    FM = np.zeros([N_tot_vars,N_tot_vars])
    for var2 in range(var1,N_tot_vars):
        start = time.clock()
        FM[var1,var2]=fisher_matrix_element(var1,var2,check_K)
        stop = time.clock()
        FM[var2,var1]=FM[var1,var2]
        np.savetxt("OUTPUT/FISHER_MATRIX_%d_K%d.csv" %(var1,CHECK_DER_K), FM)
        print "(%d, %d) --> %g (%g min.)" %(var1,var2,FM[var1,var2],(stop-start)/60.)
    return

# include "libraries/plot_stuff.pxi"

# ###########################
# Plot k-derivative comparison:
# ###########################
def plot_k_der_compare(int var_num, bin1, bin2, double k_min=0.001, double k_max=0.5, int N_k=1000):
    # N_k, N_mu = 1000, 10
    k_vect = np.linspace(k_min,k_max,N_k)
    samples = np.empty(N_k)
    samples2 = np.empty(N_k)
    # F_k.function = &argument_k
    # F_mu.function = &argument_mu
    # sdfsdf[0] = 3.
    mu = 1.
    for i_k in range(N_k):
        samples[i_k] = windowed_DER_k(k_vect[i_k],mu,bin1,bin2,var_num)
        samples2[i_k] = zero_spectrum_der_k(k_vect[i_k])*k_der(mu, k_vect[i_k], bin1, var_num)
    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(k_vect,samples,'r-',label="Windowed")
    ax1.plot(k_vect,samples2,'b-',label="Normal")
    ax1.grid(True)
    ax1.legend(loc='best')
    # ax1.set_yscale('log')
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
    fig2.savefig('plots/compare_derivatives_bins%d%d.pdf' %(bin1,bin2))
    # np.savetxt("output/FM_arg_k_integral_vars%d%d.csv" %(var1,var2),np.column_stack((k_vect,samples)))
    return

