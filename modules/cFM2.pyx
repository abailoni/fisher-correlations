###########
# Brunch:
###########
# - Trace is computed totally and only one fixed k_max is used for all bins



import sys
import numpy as np
import sympy as sym
import scipy as scy
import matplotlib.pyplot as pl
import os
import time

from os import walk
from re import match
PATH_PLOTS = "plots/"

# Cython:
cimport cython
from GSL_library cimport *
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double m)
    double sin(double x)
    double cos(double x)

# OWN LIBRARIES:
sys.path.insert(0, './modules/libraries')
from sympy_extensions import * # Some extensions to sympy library


#*************************************
# Wrappers of GSL routines in Cython:
#*************************************
include "libraries/GSL_functions.pxi"

#****************************************************************
# Import spectra from CLASS, normalise and take num. derivatives
#****************************************************************
ref_values = {'Omega_m': 0.25, 'h': 0.7, 'Omega_b': 0.0445, 'n_s': 0.967}
include "libraries/import_CLASS.pxi"


#**************************************************************
# Importing windowed data and interpolation: (with GSL)
#**************************************************************

# Tools for GLS interpolation:
cdef:
    interpolation_tools integral_1[N_bins][N_bins][N_vars]
    interpolation_tools integral_2[N_bins][N_bins][N_vars]
    interpolation_tools integral_3[N_bins][N_bins][N_vars]

# Read file names:
PATH_INPUT_WIND_SPECTRA = "INPUT/modified_spectra/"
for (dirpath, dirnames, filenames) in walk(PATH_INPUT_WIND_SPECTRA):
    names = [ fi for fi in filenames if fi.endswith(".csv") ]
    break
# Importing files:
print "\nImporting and interpolating windowed spectra. . .",
j, count1, count23 = 0, 0, 0
for filename in names:
    # Some log:
    if j%10==0:
        print".",
    j+=1
    # Determine which kind of data there is in the file:
    m = match(r"([0-9])([0-9])-var([0-9])[.\_]([a-z])", filename)
    bin1, bin2, var_num = int(m.group(1)), int(m.group(2)), int(m.group(3))
    check = m.group(4)
    # Interpolate: (var_num from )
    file_data = np.loadtxt(open(PATH_INPUT_WIND_SPECTRA+filename,"rb"),delimiter=" ")
    k_data = np.array(file_data[:,0])
    if (check=="i"):
        data_2, data_3 = np.array(file_data[:,1]), np.array(file_data[:,2])
        alloc_interp_GSL( k_data, data_2, &integral_2[bin1][bin2][var_num])
        alloc_interp_GSL( k_data, data_3, &integral_3[bin1][bin2][var_num])
        count23+=1
    else:
        data_y = np.array(file_data[:,1])
        alloc_interp_GSL( k_data, data_y, &integral_1[bin1][bin2][var_num])
        count1+=1
    # Save max_k and min_k:
    if j==1:
        k_min, k_max = k_data[0], k_data[-1]
N_expected = N_bins*(N_bins+1)/2*(N_vars+1) # For each var4 we have two files
if ( count1+count23!=N_expected ):
    print "\n\n !!!!  ---> Some windowed spectra missing: int1=%d, int23=%d --> Expected: (%d, %d)  <---  !!!!" %(count1,count23,N_bins*(N_bins+1)/2*5,N_bins*(N_bins+1)/2*3)
else:
    print "Done!"



#*************************************************
# Analytical background and derivatives: (sympy)
# compute and store data for EUCLID
#*************************************************
include "libraries/analytical_fcts.pxi"


# Compute k_max at each z: (actually here x_max is already k_max...)
print " - k_max at each z. . .",
def sigma_8(x_max, zx):
    integral_k = quad(lambda kp: kp**2 * zero_spectrum(kp) * Fourier_W_k(kp*np.pi/(2*x_max))**2, k_min, k_max, epsrel=1e-3)[0]
    integral_z = np.exp( -2. * NIntegrate( Omega_m_z**gamma/(1+z), z, 0., zx , [ref_values['gamma'], ref_values['w_1'], ref_values['w_p'], ref_values['Omega_m']] ))
    return 1/(2.*np.pi**2) * integral_z * integral_k

def sigma_8_eq(x_max, zx):
    return np.sqrt(sigma_8(x_max, zx))-np.sqrt(0.35)

from scipy.optimize import newton
def root_k_max(zx):
    return newton(sigma_8_eq,np.pi/40.,args=tuple([zx]))

cdef double[::1] k_max_data = np.array([root_k_max(zx) for zx in z_avg])
print " --> Ok\n"


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

cdef double windowed_DER(double k, int bin1, int bin2, int var):
    if bin1<=bin2:
        return eval_interp_GSL(k, &integral_1[bin1][bin2][var+1])
    else:
        return eval_interp_GSL(k, &integral_1[bin2][bin1][var+1])

# WRONG!!!
cdef double integral_B(double k, double mu, int bin1, int bin2, int var):
    return ( delete_volume_term * missing_factor*(mu*mu*eval_interp_GSL(k, &integral_2[bin1][bin2][var+1]) + eval_interp_GSL(k, &integral_3[bin1][bin2][var+1]) )/windowed_Spectrum(k,bin1,bin2) )


#--------------------------------------------------------------
# Contructing the 16 final derivatives for the Fisher Matrix:
#--------------------------------------------------------------
# Observed spectrum: (optimized!)
cdef double observed_spectrum(int bin1, int bin2, double k, double mu):
    return( Growth_bins[bin1]*Growth_bins[bin2] * bias_bins[bin1]*bias_bins[bin2]* (1+beta_bins[bin1]*mu**2)*(1+beta_bins[bin2]*mu**2) * windowed_Spectrum(k,bin1,bin2) )


# h, Omega_b and n_s: (optimized!)
cdef double der_type_A(int bin1, int bin2, double k,double mu, int var_num):
    return  observed_spectrum(bin1, bin2, k, mu)*windowed_DER(k, bin1, bin2, var_num)/windowed_Spectrum(k,bin1,bin2)

# Omega_m, w_p and w_1: (optimized!)
cdef double der_type_B(int bin1, int bin2, double k, double mu, int var_num):
    cdef double CLASS_term
    if var_num==3: #Omega_m
        CLASS_term = windowed_DER(k, bin1, bin2, var_num)/windowed_Spectrum(k,bin1,bin2)
    else:
        CLASS_term = 0
    cdef np.intp_t avg_bin = (bin1+bin2)/2
    cdef double beta_term = 1./(1+beta_bins[bin1]*mu**2)*( 2*beta_bins[bin1]*mu*mu_der(mu,avg_bin,var_num) + mu**2*Beta_der_data[var_num][bin1] ) +  1./(1+beta_bins[bin2]*mu**2)*( 2*beta_bins[bin2]*mu*mu_der(mu,avg_bin,var_num)+ mu**2*Beta_der_data[var_num][bin2] )
    # Pay attention to lnH_der_data that are computed in z_avg....
    cdef double temp = observed_spectrum(bin1, bin2, k, mu) * ( CLASS_term + lnG_der_data[var_num][bin1]+lnG_der_data[var_num][bin2] + lnH_der_data[var_num][avg_bin] - 2*lnD_der_data[var_num][avg_bin] + beta_term  )
    return(temp)

# Gamma: (optimized!)
cdef double der_gamma(int bin1, int bin2, double k, double mu):
    # Pay attention to lnH_der_data that are computed in z_avg....
    return(observed_spectrum(bin1, bin2, k, mu) * (lnG_der_data[6][bin1]+lnG_der_data[6][bin2] + 1./(1+beta_bins[bin1]*mu**2)*(mu**2*Beta_der_data[6][bin1]) + 1./(1+beta_bins[bin2]*mu**2)*(mu**2*Beta_der_data[6][bin2])) )

# Bias: (9 derivatives) (bad optimized....)
cdef double der_bias(int bin1, int bin2, double k, double mu, int bin_bias) except -1:
    bias_term = lambda bin: 1/bias_bins[bin] - 1./(1+beta_bins[bin]*mu**2)*mu**2 * Omega_m_z_fct(z_avg[bin], ref_values['w_1'], ref_values['w_p'], ref_values['Omega_m'])**ref_values['gamma'] /(bias_bins[bin]**2)
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

# Defining variables:
P_der_1, P_der_2 = np.zeros([N_bins, N_bins]), np.zeros([N_bins, N_bins])
C = np.zeros([N_bins, N_bins])
cdef:
    double[:,::1] inverse_C_v = np.empty([N_bins, N_bins])
    double[:,::1] C_v = C
    double[:,::1] P_der_1_v = P_der_1
    double[:,::1] P_der_2_v = P_der_2

# Noise matrix N:
efficiency = 0.5 # Boh...
n_dens = efficiency*np.array([3.56, 2.42, 1.81, 1.44, 0.99, 0.55, 0.29, 0.15])*1e-3
inv_n_dens = 1./n_dens
N = np.identity(N_bins) * inv_n_dens
cdef double[:,::1] N_v = N

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
cdef double argument_mu(double mu, void *input): #k, var1, var2
    cdef:
        double *params = <double*> input
        double k = params[0]
    return( k*k * trace(k,mu,<int>params[1],<int>params[2]) )

cdef double argument_k(double k, void *input): #var1, var2
    cdef double params[3]
    cdef double *vars = <double*>input
    cdef gsl_function F_mu
    F_mu.function = &argument_mu
    params[0], params[1], params[2] = k, vars[0], vars[1]
    cdef double result = eval_integration_GSL(0., 1., abs_prec, rel_prec, params, W_mu, &F_mu, MAX_ALLOC)
    # print "(%g,   %g) "%(k, result)
    return 2*result

#------------------------
# FISHER MATRIX element:
#------------------------
# new_k_min = 0.010084404404404404 #Cut the first spike
# new_k_min = 0.014284404404404404 #Cut all the first huge incomplete wave
new_k_min = 0.001
new_k_max = 0.198
def fisher_matrix_element(int var1, int var2):
    cdef gsl_function F_k
    F_k.function = &argument_k
    # Let's laugh... :/ --> IT WORKS!! :D (Apparently)
    cdef double params[2]
    params[0], params[1] = var1, var2
    cdef double result = 1./(8*np.pi**2) * eval_integration_GSL(new_k_min, new_k_max, abs_prec, rel_prec, params, W_k, &F_k, MAX_ALLOC_K)
    return(result)


#------------------------
# Computation of FM:
#------------------------
N_tot_vars = 15
def FM():
    FM = np.zeros([N_tot_vars,N_tot_vars])
    for var1 in range(N_tot_vars):
        for var2 in range(var1,N_tot_vars):
            start = time.clock()
            FM[var1,var2]=fisher_matrix_element(var1,var2)
            stop = time.clock()
            FM[var2,var1]=FM[var1,var2]
            np.savetxt("OUTPUT/FM_corr.csv", FM)
            print "(%d, %d) --> %g (%g sec.)" %(var1,var2,FM[var1,var2],(stop-start))
    return FM

def FM_var(int var1):
    FM = np.zeros([N_tot_vars,N_tot_vars])
    for var2 in range(var1,N_tot_vars):
        start = time.clock()
        FM[var1,var2]=fisher_matrix_element(var1,var2)
        stop = time.clock()
        FM[var2,var1]=FM[var1,var2]
        np.savetxt("OUTPUT/FISHER_MATRIX_%d.csv" %(var1), FM)
        print "(%d, %d) --> %g (%g min.)" %(var1,var2,FM[var1,var2],(stop-start)/60.)
    return

include "libraries/plot_stuff.pxi"

