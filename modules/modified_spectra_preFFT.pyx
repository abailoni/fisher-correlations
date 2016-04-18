

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
    # Integral DER:
    double rel_prec_kp_DER = 1e-3
    double abs_prec_kp_DER = 1e-8
    double rel_prec_z_DER = 1e-4
    double abs_prec_z_DER = 1e-10
    # Others:
    # double kp_max = 0.56, kp_min = 2e-4
    double delta_k_int_1 = 1e-1
    double delta_kp_intDER = 1.5e-1
cdef:
    gsl_integration_cquad_workspace *W_kp
    gsl_integration_workspace * W_z
W_kp = gsl_integration_cquad_workspace_alloc(MAX_ALLOC)
W_z = gsl_integration_workspace_alloc(MAX_ALLOC)



# Arguments of the integrals:
gsl_set_error_handler_off()

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


def integral_1_DEF(int num_var, int bin1, int bin2, vect_k = np.concatenate((np.linspace(k_min,0.2,130), np.linspace(0.2,0.5,60)[1:]))): # num_var [0-4]

    # Integration variables:
    cdef:
        double params_int[4]
        double extrem_1, extrem_2
        gsl_function F_kp
    params_int[1], params_int[2], params_int[3] =  bin1, bin2, num_var
    F_kp.function = &arg_integral_1

    integral = np.zeros(vect_k.shape[0])
    zero_spectr = np.zeros(vect_k.shape[0])

    cdef:
        double k
        size_t n_eval
        int count=0
    for i in range(vect_k.shape[0]):
        # Restrict extremes integral for i!=j:
        extrem_1, extrem_2 = k_min, k_max
        k = vect_k[i]
        if (bin1!=bin2):
            if (k-delta_k_int_1 > k_min ):
                extrem_1 = k- delta_k_int_1
            if (k+delta_k_int_1 < k_max ):
                extrem_2 = k + delta_k_int_1
        # Integration:
        params_int[0] = k
        if num_var==0:
            zero_spectr[i] = zero_spectrum(k)
        else:
            zero_spectr[i] = CLASS_der(k, num_var)
        if i>5:
            average = np.average( integral[i-5:i] )
        else:
            average = 42.
        if (k>0.07 and abs(average)<0.1):
            integral[i]=0
        else:
            start = time.clock()
            integral[i]  = sqrt(vol_shell_original(bin1)*vol_shell_original(bin2)) / (4*np.pi**2) * eval_integ_GSL_cquad(extrem_1, extrem_2, abs_prec_kp, rel_prec_kp, params_int, W_kp, &F_kp, &n_eval)
            stop = time.clock()
            print "(%d,%d) --> (%g, %g, %g) --> %g sec." %(bin1,bin2,k,integral[i],zero_spectr[i],stop-start)
        if (count%10==0):
                np.savetxt(PATH_WINDOWED_SPECTRA+"%d_%d-var%d_int1.csv" %(bin1,bin2,num_var),np.column_stack((vect_k,integral,zero_spectr)))
        count+=1
    np.savetxt(PATH_WINDOWED_SPECTRA+"%d_%d-var%d_int1.csv" %(bin1,bin2,num_var),np.column_stack((vect_k,integral,zero_spectr)))
    return integral

def compute_int1(var, bin1, bin2):
    # In order to get a smooth sampling of the function, the sampling
    # points k depends on bin1 and bin2.
    bin_diff = abs(bin1-bin2) # from 0 to N_bins-1
    if bin1==bin2 or bin_diff==1:
        vect_k = np.concatenate((np.linspace(k_min,0.2,130), np.linspace(0.2,0.5,60)[1:]))
    else:
        max1, max2 = 0.01, 0.1
        x = (bin_diff-2.)/(N_bins-3.) # 0 for similar bins, 1 for distant ones
        max1 = max1*(4-3*x) # max1 for distant bins, 4*max1 for similar ones
        max2 = max2*(4-3*x)
        vect_k = np.concatenate((np.linspace(k_min,max1,60), np.linspace(max1,max2,60)[1:], np.linspace(max2,0.5,20)[1:]))
    integral1 = integral_1_DEF(var,bin1,bin2,vect_k)

def wrapper_int1(var, bin1, bin2,kmin,kmax,Nk):
    vect_k=np.linspace(kmin,kmax,Nk)
    integral = integral_1_DEF(var,bin1,bin2,vect_k)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(vect_k,integral,'r-',label="Int_DER")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/int1/%d%d_var%d.pdf' %(bin1,bin2,var))
    return


# --------------------------------------------------------
# TERM WITH THE DERIVATIVE OF K_ij:
# --------------------------------------------------------

cdef double arg_angle_int_derK(double z, void *inputs): # k, kp, bin1, bin2
    cdef :
        double *params = <double*> inputs
        double k = params[0], kp = params[1]
        int bin1 = <int> params[2], bin2 = <int> params[3]
    return W2_x_derW1(k,kp,bin1,bin2,z)

cdef double arg_int_derK(double kp, void *inputs): #k, bin1, bin2
    cdef:
        double *params = <double*> inputs
        double k = params[0], params_int[4]
        int bin1 = <int> params[1], bin2 = <int> params[2]
    # Integration:
    params_int[0], params_int[1], params_int[2], params_int[3] = k, kp, bin1, bin2
    cdef gsl_function F_z
    F_z.function = &arg_angle_int_derK
    cdef double integral
    # print "  (%g  %g)" %(kp,k-kp),
    integral = eval_integration_GSL(-1., 1., abs_prec_z_DER,rel_prec_z_DER, params_int, W_z, &F_z, MAX_ALLOC)
    # print "-->%g   " %integral
    # print k
    return zero_spectrum(kp) * kp*kp * integral


def integral_DER_K_DEF(int bin1, int bin2, vect_k=np.concatenate((np.linspace(k_min,0.1,80), np.linspace(0.1,0.3,65)[1:], np.linspace(0.3,0.5,40)[1:]))):
    # Integration variables:
    cdef:
        double params_int[3]
        double extrem_1, extrem_2
        gsl_function F_kp
    params_int[1], params_int[2] =  bin1, bin2
    F_kp.function = &arg_int_derK

    # vect_k = np.concatenate((np.linspace(k_min,0.1,80), np.linspace(0.1,0.3,65)[1:], np.linspace(0.3,0.5,40)[1:]))
    integral = np.zeros(vect_k.shape[0])

    cdef:
        double k
        size_t n_eval
        int count=0
    for i in range(vect_k.shape[0]):
        # Restrict extremes integral for i!=j:
        extrem_1, extrem_2 = k_min, k_max
        if (bin1!=bin2):
            if (k-delta_kp_intDER > k_min ):
                extrem_1 = k- delta_kp_intDER
            if (k+delta_kp_intDER < k_max ):
                extrem_2 = k + delta_kp_intDER

        k = vect_k[i]
        # Integration:
        params_int[0] = k
        start = time.clock()
        integral[i]  = eval_integ_GSL_cquad(extrem_1, extrem_2, abs_prec_kp_DER, rel_prec_kp_DER, params_int, W_kp, &F_kp, &n_eval)
        stop = time.clock()
        print "BINS: (%d,%d) || (k=%g, int_DER=%g) --> %g sec." %(bin1,bin2,k,integral[i],stop-start)
        if (count%10==0):
            np.savetxt(PATH_WINDOWED_SPECTRA+"%d_%d_intDER.csv" %(bin1,bin2),np.column_stack((vect_k,integral)))
        count+=1
    np.savetxt(PATH_WINDOWED_SPECTRA+"%d_%d_intDER.csv" %(bin1,bin2),np.column_stack((vect_k,integral)))
    return integral


def compute_intDER(bin1,bin2):
    # In order to get a smooth sampling of the function, the sampling
    # points k depends on bin1 and bin2.
    global abs_prec_z_DER
    if (bin1!=bin2):
        abs_prec_z_DER = 1e-12
        k_distant = np.concatenate((np.linspace(k_min,0.004,100), np.linspace(0.004,0.01,70)[1:], np.linspace(0.01,0.05,200)[1:]))
        k_close = k_distant*4 # bin1 similar to bin2
        bin_diff = abs(bin1-bin2) # from 1 to N_bins-1
        x = (bin_diff-1.)/(N_bins-2.) # 0 for similar bins, 1 for distant ones
        vect_k = k_distant + (k_close-k_distant)*(1-x)
        # Correction for close bins:
        if (bin_diff==1):
            vect_k = np.concatenate(( vect_k, np.linspace(0.2,0.3,40)[1:]))
    else:
        abs_prec_z_DER = 1e-10
        Nk_1, Nk_2, Nk_3, Nk_4, Nk_5, Nk_6 = 30, 60, 80, 50, 20, 10
        if (bin1 > N_bins/3): # last bins are slower
            Nk_2, Nk_3, Nk_4, Nk_5 = Nk_2/2, Nk_3/3, Nk_4/2, Nk_5/3*2
        vect_k = np.concatenate((np.linspace(k_min,0.02,Nk_1), np.linspace(0.02,0.1,Nk_2)[1:], np.linspace(0.1,0.2,Nk_3)[1:], np.linspace(0.2,0.3,Nk_4)[1:], np.linspace(0.3,0.4,Nk_5)[1:], np.linspace(0.4,0.5,Nk_6)[1:] ))
    integral = integral_DER_K_DEF(bin1,bin2,vect_k)






def wrapper(bin1,bin2,kmin,kmax,Nk):
    vect_k=np.linspace(kmin,kmax,Nk)
    integral = integral_DER_K_DEF(bin1,bin2,vect_k)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(vect_k,integral,'r-',label="Int_DER")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/intDER/%d%d.pdf' %(bin1,bin2))
    return



def plot_angDER(bin1,bin2,k,kp, N_z = 600,min_z=-1, max_z=1):
    z_vect = np.linspace(min_z,max_z,N_z)
    samples = np.zeros(N_z)
    cdef double params[4]
    params[0], params[1], params[2], params[3]= k, kp, bin1, bin2
    for i in range(N_z):
            samples[i] = arg_angle_int_derK(z_vect[i],params)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(z_vect,samples,'r-',label="Arg int. z --> DER K_ij")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/intDER/arg_angle_intDER_%d%d_%g_%g.pdf' %(bin1,bin2,k,kp))
    return

def plot_intDER(bin1,bin2, k, N_kp = 600, kp_start=k_min, kp_stop=k_max, rel_input=rel_prec_z_DER, abs_input=abs_prec_z_DER):
    global abs_prec_z_DER, rel_prec_z_DER
    abs_prec_z_DER = abs_input
    rel_prec_z_DER = rel_input
    kp_vect = np.linspace(kp_start,kp_stop,N_kp)
    samples = np.zeros(N_kp)
    cdef double params[3]
    params[0], params[1], params[2]= k, bin1, bin2
    start = time.clock()
    for i in range(N_kp):
            samples[i] = arg_int_derK(kp_vect[i],params)
    print time.clock()-start

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(kp_vect,samples,'r-',label="Arg. int. kp --> DER K_ij")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/intDER/arg_intDER_%d%d_%g.pdf' %(bin1,bin2,k))
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


def plot_int1(bin1,bin2, k, num_var, N_kp = 600, kp_start=k_min, kp_stop=k_max):
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

# def plot_der_W1(bin1,bin2,k,kp, num_var, N_z = 600, min_z=-1, max_z=1):
#     z_vect = np.linspace(min_z,max_z,N_z)
#     samples = np.zeros(N_z)
#     for i in range(N_z):
#             samples[i] = der_W1(k,kp,bin1,bin2,z_vect[i])
#     fig2=pl.figure()
#     ax1=fig2.add_subplot(111)
#     ax1.plot(z_vect,samples,'r-',label="der_W1")
#     # ax1.set_xlabel("$k$ [$h$/Mpc]")
#     # ax1.set_ylabel("[Mpc/$h$]^3")
#     ax1.grid(True)
#     ax1.legend(loc='best')
#     fig2.savefig('plots/der_W1_%d%d_var%d_%g_%g.pdf' %(bin1,bin2,num_var,k,kp))
#     return

# def plot_der_Wcompiled(bin, N_kmod = 1000, start=0, end=0.5):
#     k_vect = np.linspace(start,end,N_kmod)
#     samples = np.zeros(N_kmod)
#     for i in range(N_kmod):
#             samples[i] = derW_compiled(k_vect[i],bin)

#     fig2=pl.figure()
#     ax1=fig2.add_subplot(111)
#     ax1.plot(k_vect,samples,'r-',label="der_W_modK")
#     # ax1.set_xlabel("$k$ [$h$/Mpc]")
#     # ax1.set_ylabel("[Mpc/$h$]^3")
#     ax1.grid(True)
#     ax1.legend(loc='best')
#     fig2.savefig('plots/der_W_modK.pdf' )
#     return

# def plot_Wcompiled(bin, N_kmod = 1000, start=0, end=0.5):
#     k_vect = np.linspace(start,end,N_kmod)
#     samples = np.zeros(N_kmod)
#     for i in range(N_kmod):
#             samples[i] = W_compiled(k_vect[i],bin)

#     fig2=pl.figure()
#     ax1=fig2.add_subplot(111)
#     ax1.plot(k_vect,samples,'r-',label="der_W_modK")
#     # ax1.set_xlabel("$k$ [$h$/Mpc]")
#     # ax1.set_ylabel("[Mpc/$h$]^3")
#     ax1.grid(True)
#     ax1.legend(loc='best')
#     fig2.savefig('plots/W_modK.pdf' )
#     return

