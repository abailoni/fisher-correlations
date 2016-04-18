# Old branch: 14-03-2016
#
# int2 and int3 are still present


# Importing modules and defining default data:

include "libraries/header.pxi"

include "libraries/import_CLASS.pxi"

include "libraries/analytical_fcts.pxi"


# Not necessary:
def all():
    import_CLASS_data()
    compute_survey_DATA()

#**************************************************************
#**************************************************************
#
# COMPUTING AND SAVING PRESENT WINDOWED SPECTRA:
#
#**************************************************************
#**************************************************************

PATH_WINDOWED_SPECTRA = "INPUT/modified_spectra/"
PATH_WINDOWED_SPECTRA_PLOTS = "plots/modified_spectra/"

# Integration variables:

cdef enum:
    max_alloc_const = 5000
cdef:
    size_t MAX_ALLOC = max_alloc_const
    # Integral 1:
    double rel_prec_kp = 1e-3
    double abs_prec_kp = 1e-8
    double rel_prec_z = 1e-5
    double abs_prec_z = 1e-14
    # Integral 2 and 3:
    double rel_prec_kp_23 = 1e-3
    double rel_prec_z_23 = 1e-4
    double abs_prec_kp_2 = 1e-8
    double abs_prec_kp_3 = 1e-8
    double abs_prec_z_2 = 1e-10
    double abs_prec_z_3 = 1e-9
    # Others:
    double safe_kp_max_int = 0.56, safe_kp_min_int = 2e-4
    double delta_k_int_1 = 8e-2
    # double delta_k_int_23 = 1e-2
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


def integral_1_DEF(int num_var, int bin1, int bin2): # num_var [0-4]
    new_k_min=k_min
    new_k_max=0.5

    # Integration variables:
    cdef:
        double params_int[4]
        double extrem_1, extrem_2
        gsl_function F_kp
    params_int[1], params_int[2], params_int[3] =  bin1, bin2, num_var
    F_kp.function = &arg_integral_1

    vect_k = np.concatenate((np.linspace(new_k_min,0.2,130), np.linspace(0.2,0.5,60)[1:]))
    integral = np.empty(vect_k.shape[0])
    zero_spectr = np.empty(vect_k.shape[0])

    cdef:
        double k
        size_t n_eval
    for i in range(vect_k.shape[0]):
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
        if (k>0.07 and abs(integral[i-1])<0.1):
            integral[i]=0
        else:
            start = time.clock()
            integral[i]  = np.sqrt(vol_shell_original(bin1)*vol_shell_original(bin2)) / (4*np.pi**2) * eval_integ_GSL_cquad(extrem_1, extrem_2, abs_prec_kp, rel_prec_kp, params_int, W_kp, &F_kp, &n_eval)
            stop = time.clock()
            print "(%d,%d) --> (%g, %g, %g) --> %g sec." %(bin1,bin2,k,integral[i],zero_spectr[i],stop-start)
    np.savetxt(PATH_WINDOWED_SPECTRA+"%d%d-var%d_int1.csv" %(bin1,bin2,num_var),np.column_stack((vect_k,integral,zero_spectr)))
    return


# --------------------------------------------------------
# TERM WITH THE DERIVATIVE OF K_ij:
# --------------------------------------------------------

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
    # if ( abs(k-kp)<1e-1):
    integral = eval_integration_GSL(-1., 1., abs_prec_z_2,rel_prec_z_23, params_int, W_z, &F_z, MAX_ALLOC)
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

    integral = eval_integration_GSL(-1., 1., abs_prec_z_3, rel_prec_z_23, params_int, W_z, &F_z, MAX_ALLOC)
    return zero_spectrum(kp) * kp*kp * integral

def integral_DER_K_DEF(int num_var, int bin1, int bin2): # num_var [4-6]
    global abs_prec_z_2, abs_prec_z_3, abs_prec_kp_2, abs_prec_kp_3

    # Integration variables:
    cdef:
        double params_int[4]
        double extrem_1, extrem_2
        gsl_function F_kp_2, F_kp_3
    params_int[1], params_int[2], params_int[3] =  bin1, bin2, num_var
    F_kp_2.function = &arg_integral_2
    F_kp_3.function = &arg_integral_3

    vect_k = np.concatenate((np.linspace(k_min,0.1,80), np.linspace(0.1,0.3,65)[1:], np.linspace(0.3,0.5,40)[1:]))
    integral2, integral3 = np.zeros(vect_k.shape[0]), np.zeros(vect_k.shape[0])

    cdef:
        double k
        size_t n_eval
        int count=0
    for i in range(vect_k.shape[0]):
        # Restrict extremes integral for i!=j:
        extrem_1, extrem_2 = safe_kp_min_int, safe_kp_max_int
        k = vect_k[i]
        # Check if increase integral precision:
        # if (k>0.1):
        #     abs_prec_kp_2=1e-8
        #     abs_prec_z_2=1e-10
        #     abs_prec_kp_3=1e-7
        if (k>0.35):
            abs_prec_kp_2=1e-9
        # Integration:
        params_int[0] = k
        start = time.clock()
        integral2[i]  = np.sqrt(vol_shell_original(bin1)*vol_shell_original(bin2)) / (4*np.pi**2) * eval_integ_GSL_cquad(extrem_1, extrem_2, abs_prec_kp_2, rel_prec_kp_23, params_int, W_kp2, &F_kp_2, &n_eval)
        stop = time.clock()
        print "BINS: (%d,%d) || (k=%g, 2=%g) --> %g sec." %(bin1,bin2,k,integral2[i],stop-start),
        start = time.clock()
        integral3[i]  = np.sqrt(vol_shell_original(bin1)*vol_shell_original(bin2)) / (4*np.pi**2) * eval_integ_GSL_cquad(extrem_1, extrem_2, abs_prec_kp_3, rel_prec_kp_23, params_int, W_kp3, &F_kp_3, &n_eval)
        stop = time.clock()
        print "|| (3=%g) --> %g sec." %(integral3[i],stop-start)
        if (count%10==0):
                np.savetxt(PATH_WINDOWED_SPECTRA+"%d%d-var%d_int23.csv" %(bin1,bin2,num_var),np.column_stack((vect_k,integral2,integral3)))
        count+=1
    np.savetxt(PATH_WINDOWED_SPECTRA+"%d%d-var%d_int23.csv" %(bin1,bin2,num_var),np.column_stack((vect_k,integral2,integral3)))
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

def plot_der_Wcompiled(bin, N_kmod = 1000, start=0, end=0.5):
    k_vect = np.linspace(start,end,N_kmod)
    samples = np.zeros(N_kmod)
    for i in range(N_kmod):
            samples[i] = derW_compiled(k_vect[i],bin)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(k_vect,samples,'r-',label="der_W_modK")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/der_W_modK.pdf' )
    return

def plot_Wcompiled(bin, N_kmod = 1000, start=0, end=0.5):
    k_vect = np.linspace(start,end,N_kmod)
    samples = np.zeros(N_kmod)
    for i in range(N_kmod):
            samples[i] = W_compiled(k_vect[i],bin)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(k_vect,samples,'r-',label="der_W_modK")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/W_modK.pdf' )
    return

