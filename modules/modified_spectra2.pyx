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
INPUT_DATA = "inputs/"


# Cython:
cimport cython
from GSL_library cimport *
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double m)
    double sin(double x)
    double cos(double x)


# ctypedef struct gsl_function:
#     double (* function) (double x, void * params)
#     void * params

# OWN LIBRARIES:
sys.path.insert(0, './libraries')
from sympy_extensions import * # Some extensions to sympy library


#*************************************
# Wrappers of GSL routines in Cython:
#*************************************
include "libraries/GSL_functions.pxi"


#**************************************************************
# Import spectra from CLASS, normalise and take num. der.:
#**************************************************************
ref_values = {'Omega_m': 0.25, 'h': 0.7, 'Omega_b': 0.0445, 'n_s': 0.967}
include "libraries/import_CLASS.pxi"


#**************************************************************
# Analytical background and derivatives: (sympy)
# compute and store data for EUCLID
#**************************************************************
include "libraries/analytical_fcts.pxi"


#**************************************************************
#**************************************************************
#
# COMPUTING AND SAVING PRESENT WINDOWED SPECTRA:
#
#**************************************************************
#**************************************************************

PATH_WINDOWED_SPECTRA = "OUTPUT/modified_spectra/"
PATH_WINDOWED_SPECTRA_PLOTS = "plots/modified_spectra/"

# Integration variables:

cdef enum:
    max_alloc_const = 5000
cdef:
    size_t MAX_ALLOC = max_alloc_const
    double rel_prec_kp = 1e-3
    double abs_prec_kp = 1e-8
    double rel_prec_z = 1e-5
    double abs_prec_z = 1e-14
    double rel_prec_kp_23 = 1e-3
    double abs_prec_kp_23 = 1e-8
    # double rel_prec_z_23 = 1e-5
    # double abs_prec_z_23 = 1e-14
    double rel_prec_z_23 = 1e-4
    double abs_prec_z_23 = 1e-5
    double safe_kp_max_int = 0.56, safe_kp_min_int = 2e-4
    double delta_k_int_1 = 8e-2
    double delta_k_int_23 = 1e-2
cdef:
    gsl_integration_cquad_workspace *W_kp, *W_kp2, *W_kp3, *W_z_prova
    gsl_integration_workspace * W_z
W_kp, W_kp2, W_kp3, W_z_prova = gsl_integration_cquad_workspace_alloc(MAX_ALLOC), gsl_integration_cquad_workspace_alloc(MAX_ALLOC), gsl_integration_cquad_workspace_alloc(MAX_ALLOC), gsl_integration_cquad_workspace_alloc(MAX_ALLOC)
W_z = gsl_integration_workspace_alloc(MAX_ALLOC)



# Arguments of the integrals:

cdef double arg_angle_integral_1(double z, void *inputs): #k, kp, bin1, bin2
    cdef:
        double *params = <double*> inputs
        double k = params[0], kp = params[1]
        int bin1 = <int> params[2], bin2 = <int> params[3]
    return( K(np.sqrt(k*k+kp*kp-2*k*kp*z),bin1,bin2) )


cdef double arg_integral_1(double kp, void *inputs): #k, bin1, bin2, num_var
    cdef:
        double *params = <double*> inputs
        double k = params[0], params_int[4]
        int bin1 = <int> params[1], bin2 = <int> params[2], num_var = <int> params[3]

    params_int[0], params_int[1], params_int[2],params_int[3] = k, kp, bin1, bin2
    cdef gsl_function F_z
    F_z.function = &arg_angle_integral_1
    cdef double integral = eval_integration_GSL(-1., 1., abs_prec_z, rel_prec_z, params_int, W_z, &F_z, MAX_ALLOC)
    if num_var==0:
        return(zero_spectrum(kp) * kp*kp * integral)
    else:
        return(CLASS_der(kp,num_var) * kp*kp * integral)


def integral_1(int num_var, int bin1, int bin2, int N_k=200, new_k_min=k_min, new_k_max=0.5):
    # Improving time: (define variables)
    cdef:
        double last_results[20], last_eval[5]
        int current_index = 0, current_index_eval = 0, check_exit = 0

    # Integration variables:
    cdef:
        double params_int[4]
        double extrem_1, extrem_2
        gsl_function F_kp
    params_int[1], params_int[2], params_int[3] =  bin1, bin2, num_var
    F_kp.function = &arg_integral_1

    vect_k = np.linspace(new_k_min,new_k_max,N_k)
    integral = np.empty(vect_k.shape[0])
    zero_spectr = np.empty(vect_k.shape[0])

    cdef:
        double k
        size_t n_eval
    for i in range(N_k):
        if (check_exit==0):
            # Restrict extremes integral for i!=j:
            extrem_1, extrem_2 = safe_kp_min_int, safe_kp_max_int
            k = vect_k[i]
            if (bin1!=bin2):
                if (k-delta_k_int_1 > safe_kp_min_int ):
                    extrem_1 = k- delta_k_int_1
                if (k+delta_k_int_1 < safe_kp_max_int ):
                    extrem_2 = k + delta_k_int_1
            # Integration:
            params_int[0] = k
            if num_var==0:
                zero_spectr[i] = zero_spectrum(k)
            else:
                zero_spectr[i] = CLASS_der(k, num_var)
            check = 0
            # rel_prec_kp_temp = rel_prec_kp
            # while (check==0):
            start = time.clock()
            integral[i]  = np.sqrt(vol_shell_original(bin1)*vol_shell_original(bin2)) / (4*np.pi**2) * eval_integ_GSL_cquad(extrem_1, extrem_2, abs_prec_kp, rel_prec_kp, params_int, W_kp, &F_kp, &n_eval)
                # if (bin1==bin2):
                #     if ( abs((zero_spectr[i]-integral[i])/zero_spectr[i])>1 and k>0.05 ):
                #         rel_prec_kp_temp= rel_prec_kp_temp/10.
                #         print rel_prec_kp_temp
                #     else:
                #         check=1
                # else:
                #     check=1
            stop = time.clock()
            print "(%d,%d) --> (%g, %g, %g) --> %g sec." %(bin1,bin2,k,integral[i],zero_spectr[i],stop-start)
    np.savetxt(PATH_WINDOWED_SPECTRA+"%d%d-var%d.csv" %(bin1,bin2,num_var),np.column_stack((vect_k,integral,zero_spectr)))
    return

def integral_A(bin1,bin2,int N_k=120):
    print "(((((%d,%d))))))" %(bin1,bin2)
    for var in range(5):
        integral_1(var,bin1,bin2,N_k)


# Integral 1 at higher momenta:
def integral_B(bin1,bin2):
    for var in range(5):
        integral_1(var,bin1,bin2,60,0.2,0.5)

# TERM WITH THE DERIVATIVE OF K_ij:

cdef double arg_angle_integral_2(double z, void *inputs): # k, kp, bin1, bin2, num_var --> var=[4-6]
    cdef :
        double *params = <double*> inputs
        double k = params[0], kp = params[1]
        int bin1 = <int> params[2], bin2 = <int> params[3], num_var = <int> params[4]
        double result = k *  ( der_W1(k, kp, bin1, bin2, z) * (lnH_der_data[num_var-1][bin1] + lnD_der_data[num_var-1][bin1]) + der_W2(k, kp, bin1, bin2, z)*(lnH_der_data[num_var-1][bin2] + lnD_der_data[num_var-1][bin2]) )
    # print "Ang=%g, z=%g, kp=%g, k=%g " %(result, z,kp,k)
    return  result

cdef:
    double[::1] z_max=np.zeros(2), z_min=np.zeros(2)
    # double *prova

cdef double arg_integral_2(double kp, void *inputs): #k, bin1, bin2, num_var
    cdef:
        double *params = <double*> inputs
        double k = params[0], params_int[5]
        int bin1 = <int> params[1], bin2 = <int> params[2], num_var = <int> params[3]
    # Integration:
    params_int[0], params_int[1], params_int[2], params_int[3], params_int[4] =k, kp, bin1, bin2, num_var
    cdef gsl_function F_z
    F_z.function = &arg_angle_integral_2
    cdef double singol[4]
    singol[0]=1.
    cdef double integral
    # if (False):
    #     integral = eval_integration_GSL(z_min[0], 0.9999999, abs_prec_z_23, rel_prec_z_23, params_int, W_z, &F_z, MAX_ALLOC)
    #     # Increasing this precision gives an error:
    #     # integral = integral + eval_integration_GSL(0.999, 1., 1e-4, 1e-3, params_int, W_z, &F_z, MAX_ALLOC)
    # else:
    integral = eval_integration_GSL(z_min[0], z_max[0], abs_prec_z_23, rel_prec_z_23, params_int, W_z, &F_z, MAX_ALLOC)
    # print "### Int=%g, kp=%g, k=%g " %(integral, kp,k)
    return zero_spectrum(kp) * kp*kp * integral



cdef double arg_angle_integral_3(double z, void *inputs): # k, kp, bin1, bin2, num_var--> var=[4-6]
    cdef :
        double *params = <double*> inputs
        double k = params[0], kp = params[1]
        int bin1 = <int> params[2], bin2 = <int> params[3], num_var = <int> params[4]
    return  - k *  ( der_W1(k, kp, bin1, bin2, z) * lnD_der_data[num_var-1][bin1] + der_W2(k, kp, bin1, bin2, z) * lnD_der_data[num_var-1][bin2] )

cdef double arg_integral_3(double kp, void *inputs): #k, bin1, bin2, num_var
    cdef:
        double *params = <double*> inputs
        double k = params[0], params_int[5]
        int bin1 = <int> params[1], bin2 = <int> params[2], num_var = <int> params[3]
    # Integration:
    params_int[0], params_int[1], params_int[2], params_int[3], params_int[4] =k, kp, bin1, bin2, num_var
    cdef gsl_function F_z
    F_z.function = &arg_angle_integral_3
    cdef double singol[4]
    singol[0]=1.
    cdef double integral
    cdef size_t n_eval_3
    # if ( abs((k-kp)/k)<1e-3):
    if (True):
        # print "A1",
        integral = eval_integration_GSL(z_min[0], 0.9999999, abs_prec_z_23, rel_prec_z_23, params_int, W_z, &F_z, MAX_ALLOC)
        # Increasing this precision gives an error:
        # print "#### ddd %g "%(abs(k-kp))
        # integral = integral + int_trapezi(0.9995,1.,1e-3,params_int,&arg_angle_integral_3,300)
        # integral = integral + eval_integration_GSL_noadp(0.9995, 1.,  1e-5, 1e-4, params_int, &F_z, &n_eval_3)
        # integral = integral + eval_integration_GSL(0.9995, 1.,  1e-5, 1e-4, params_int, W_z, &F_z, MAX_ALLOC)
        # integral = integral + eval_integ_GSL_cquad(0.9995, 1., 1e-3, 1e-2, params_int, W_z_prova, &F_z, &n_eval_3)
        # print integral
    else:
        # print "B",
        integral = eval_integration_GSL(z_min[0], z_max[0], abs_prec_z_23, rel_prec_z_23, params_int, W_z, &F_z, MAX_ALLOC)
    # cdef double integral = eval_integration_GSL(z_min[0], z_max[0], abs_prec_z_23, rel_prec_z_23, params_int, W_z, &F_z, MAX_ALLOC)
    return zero_spectrum(kp) * kp*kp * integral


def integral_DER_K(int num_var, int bin1, int bin2, int N_k=200, new_k_min = k_min, new_z_min=-1., new_z_max=1.):
    z_min[0] = new_z_min
    z_max[0] = new_z_max

    # Improving time: (define variables)
    cdef:
        double last_results[20], last_eval[5]
        int current_index = 0, current_index_eval = 0, check_exit = 0

    # Integration variables:
    cdef:
        double params_int[4]
        double extrem_1, extrem_2
        gsl_function F_kp_2, F_kp_3
    params_int[1], params_int[2], params_int[3] =  bin1, bin2, num_var
    F_kp_2.function = &arg_integral_2
    F_kp_3.function = &arg_integral_3



    vect_k = np.linspace(new_k_min,0.2,N_k)
    integral2, integral3 = np.empty(vect_k.shape[0]), np.empty(vect_k.shape[0])
    zero_spectr = np.empty(vect_k.shape[0])

    cdef:
        double k
        size_t n_eval
    for i in range(N_k):
        if (check_exit==0):
            # Restrict extremes integral for i!=j:
            extrem_1, extrem_2 = safe_kp_min_int, safe_kp_max_int
            k = vect_k[i]
            # if (bin1!=bin2):
            #     if (k-delta_k_int_23 > safe_kp_min_int ):
            #         extrem_1 = k- delta_k_int_23
            #     if (k+delta_k_int_23 < safe_kp_max_int ):
            #         extrem_2 = k + delta_k_int_23
            # Integration:
            params_int[0] = k
            zero_spectr[i] = zero_spectrum(k)
            check = 0
            # rel_prec_kp_temp = rel_prec_kp
            # while (check==0):
            start = time.clock()
            integral2[i]  = np.sqrt(vol_shell_original(bin1)*vol_shell_original(bin2)) / (4*np.pi**2) * eval_integ_GSL_cquad(extrem_1, extrem_2, abs_prec_kp_23, rel_prec_kp_23, params_int, W_kp2, &F_kp_2, &n_eval)
                # if (bin1==bin2):
                #     if ( abs((zero_spectr[i]-integral[i])/zero_spectr[i])>1 and k>0.05 ):
                #         rel_prec_kp_temp= rel_prec_kp_temp/10.
                #         print rel_prec_kp_temp
                #     else:
                #         check=1
                # else:
                #     check=1
            stop = time.clock()
            print "BINS: (%d,%d) || (k=%g, 2=%g) --> %g sec." %(bin1,bin2,k,integral2[i],stop-start),
            start = time.clock()
            # integral3[i]=0
            integral3[i]  = np.sqrt(vol_shell_original(bin1)*vol_shell_original(bin2)) / (4*np.pi**2) * eval_integ_GSL_cquad(extrem_1, extrem_2, abs_prec_kp, rel_prec_kp, params_int, W_kp3, &F_kp_3, &n_eval)
                # if (bin1==bin2):
                #     if ( abs((zero_spectr[i]-integral[i])/zero_spectr[i])>1 and k>0.05 ):
                #         rel_prec_kp_temp= rel_prec_kp_temp/10.
                #         print rel_prec_kp_temp
                #     else:
                #         check=1
                # else:
                #     check=1
            stop = time.clock()
            print "|| (3=%g) --> %g sec." %(integral3[i],stop-start)
    np.savetxt(PATH_WINDOWED_SPECTRA+"%d%d-var%d_DER_int.csv" %(bin1,bin2,num_var),np.column_stack((vect_k,integral2,integral3,zero_spectr)))
    return


def plot_ang2(bin1,bin2,k,kp, num_var, N_z = 600,min_z=-1, max_z=1):
    z_vect = np.linspace(min_z,max_z,N_z)
    samples = np.zeros(N_z)
    cdef double params[5]
    params[0], params[1], params[2], params[3], params[4]= k, kp, bin1, bin2, num_var
    for i in range(N_z):
        samples[i] = arg_angle_integral_2(z_vect[i],params)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(z_vect,samples,'r-',label="Arg int. z --> DER K_ij")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/prove/arg_angle_integral_2_%d%d_var%d_%g_%g.pdf' %(bin1,bin2,num_var,k,kp))
    return

def plot_ang3(bin1,bin2,k,kp, num_var, N_z = 600,min_z=-1, max_z=1):
    z_vect = np.linspace(min_z,max_z,N_z)
    samples = np.zeros(N_z)
    cdef double params[5]
    params[0], params[1], params[2], params[3], params[4]= k, kp, bin1, bin2, num_var
    for i in range(N_z):
        samples[i] = arg_angle_integral_3(z_vect[i],params)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(z_vect,samples,'r-',label="Arg int. z --> DER K_ij")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/prove/arg_angle_integral_3_%d%d_var%d_%g_%g.pdf' %(bin1,bin2,num_var,k,kp))
    return

def plot_ang1(bin1,bin2,k,kp, num_var,  N_z = 600,min_z=-1, max_z=1):
    z_vect = np.linspace(min_z,max_z,N_z)
    samples = np.zeros(N_z)
    cdef double params[4]
    params[0], params[1], params[2], params[3]= k, kp, bin1, bin2
    for i in range(N_z):
            samples[i] = arg_angle_integral_1(z_vect[i],params)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(z_vect,samples,'r-',label="Arg int. z --> Main term")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/prove/arg_angle_integral_1_%d%d_var%d_%g_%g.pdf' %(bin1,bin2,num_var,k,kp))
    return

def plot_int2(bin1,bin2, k, num_var, N_kp = 600, kp_start=safe_kp_min_int, kp_stop=safe_kp_max_int):
    kp_vect = np.linspace(kp_start,kp_stop,N_kp)
    samples = np.zeros(N_kp)
    cdef double params[4]
    params[0], params[1], params[2], params[3]= k, bin1, bin2, num_var
    for i in range(N_kp):
            samples[i] = arg_integral_2(kp_vect[i],params)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(kp_vect,samples,'r-',label="Arg. int. kp --> DER K_ij")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/arg_integral_2_%d%d_var%d_%g.pdf' %(bin1,bin2,num_var,k))
    return

def plot_int1(bin1,bin2, k, num_var, N_kp = 600, kp_start=safe_kp_min_int, kp_stop=safe_kp_max_int):
    kp_vect = np.linspace(kp_start,kp_stop,N_kp)
    samples = np.zeros(N_kp)
    cdef double params[4]
    params[0], params[1], params[2], params[3]= k, bin1, bin2, num_var
    for i in range(N_kp):
            samples[i] = arg_integral_1(kp_vect[i],params)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(kp_vect,samples,'r-',label="Arg. int. kp --> Main term")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/arg_integral_1_%d%d_var%d_%g.pdf' %(bin1,bin2,num_var,k))
    return

def plot_der_W1(bin1,bin2,k,kp, num_var, N_z = 600, min_z=-1, max_z=1):
    z_vect = np.linspace(min_z,max_z,N_z)
    samples = np.zeros(N_z)
    for i in range(N_z):
            samples[i] = der_W1(k,kp,bin1,bin2,z_vect[i])

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(z_vect,samples,'r-',label="der_W1")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/der_W1_%d%d_var%d_%g_%g.pdf' %(bin1,bin2,num_var,k,kp))
    return


