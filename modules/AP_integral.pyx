


# Importing modules and defining default data:

include "libraries/header.pxi"

# include "libraries/import_CLASS.pxi"

include "libraries/CAMB.pxi"

include "libraries/analytical_fcts.pxi"

# It does all the first necessary things:
def init():
    compute_CAMB_spectra()
    compute_survey_DATA()
    # store_conv_spectra()

#**************************************************************
#**************************************************************
#
# COMPUTING AND SAVING PRESENT WINDOWED SPECTRA:
#
#**************************************************************
#**************************************************************

PATH_WINDOWED_SPECTRA = "INPUT/AP_integrals/"
# PATH_WINDOWED_SPECTRA_PLOTS = "plots/modified_spectra/"

# Integration variables:

cdef enum:
    max_alloc_const = 5000
cdef:
    size_t MAX_ALLOC = max_alloc_const
    # Integral AP:
    double rel_prec_k1 = 1e-3
    double abs_prec_k1 = 1e-8
    double rel_prec_mu1 = 1e-5
    double abs_prec_mu1 = 1e-14
    # double rel_prec_phi1 = 1e-3
    # double abs_prec_phi1 = 1e-8

    # Others:
    # double kp_max = 0.56, kp_min = 2e-4
    double delta_k_int_1 = 1e-1
    double delta_kp_intDER = 1.5e-1
cdef:
    gsl_integration_cquad_workspace *W_k1
    gsl_integration_workspace * W_mu1
    gsl_integration_workspace * W_phi1
W_k1 = gsl_integration_cquad_workspace_alloc(MAX_ALLOC)
W_mu1 = gsl_integration_workspace_alloc(MAX_ALLOC)
W_phi1 = gsl_integration_workspace_alloc(MAX_ALLOC)



# Arguments of the integrals:
gsl_set_error_handler_off()

# cdef double arg_integral_phi1(double phi1, void *inputs): #k, k1, mu, mu1, bin1, _phi1bin2
#     cdef:
#         double *params = <double*> inputs
#         double k = params[0], k1 = params[1]
#         double mu = params[2], mu1 = params[3]
#         int bin1 = <int> params[4], bin2 = <int> params[5]
#     cdef double square_root = sqrt(k*k + k1*k1 + 2*k*k1* (mu*mu1 + cos(phi1)*sqrt( (1-mu*mu)*(1-mu1*mu1) )  ))
#     return K(square_root,bin1,bin2)


cdef double arg_integral_mu1(double mu1, void *inputs): #k, k1, mu, bin1, bin2, num_var
    cdef:
        double *params = <double*> inputs
        double k = params[0], k1 = params[1], mu = params[2]
        int bin1 = <int> params[3], bin2 = <int> params[4], num_var = <int> params[5]

    # This works only with bin1=bin2:
    cdef double integral_phi1 = Lambda_Ev(phi1_intregral_data, k1, mu1, mu, lnH_der_data[num_var][bin1], lnD_der_data[num_var][bin1])
    cdef double square_root = sqrt(k*k + k1*k1 -2*k*k1*mu1)
    return integral_phi1 * K(square_root,bin1,bin2)




cdef double arg_integral_k1(double k1, void *inputs): #k, mu, bin1, bin2, num_var
    cdef:
        double *params = <double*> inputs
        double k = params[0], mu = params[1]
        int bin1 = <int> params[2], bin2 = <int> params[3], num_var = <int> params[4]
        double params_int_mu1[6]

    params_int_mu1[0], params_int_mu1[1], params_int_mu1[2],params_int_mu1[3],params_int_mu1[4],params_int_mu1[5] = k, k1, mu, bin1, bin2, num_var
    cdef gsl_function F_mu1
    F_mu1.function = &arg_integral_mu1
    cdef double integral_mu1 = eval_integration_GSL(-1., 1., abs_prec_mu1, rel_prec_mu1, params_int_mu1, W_mu1, &F_mu1, MAX_ALLOC)

    return zero_spectrum_der_k(k1) * k1*k1 * integral_mu1


def AP_integral(int num_var, int bin1, int bin2, vect_k = np.concatenate((np.linspace(k_min,0.2,130), np.linspace(0.2,0.5,60)[1:])), vect_mu=np.linspace(-1.,1.,5)): # num_var [2-4]

    # Integration variables:
    cdef:
        double params_int_k1[5]
        double extrem_1, extrem_2
        gsl_function F_k1
    params_int_k1[2], params_int_k1[3], params_int_k1[4] =  bin1, bin2, num_var
    F_k1.function = &arg_integral_k1

    integral = np.zeros((vect_mu.shape[0],vect_k.shape[0]))

    cdef:
        double k, mu
        size_t n_eval
        int count=0
    for j_mu in range(vect_mu.shape[0]):
        mu = vect_mu[j_mu]
        params_int_k1[1] = mu
        for i in range(vect_k.shape[0]):
            # Restrict extremes integral for i!=j:
            extrem_1, extrem_2 = k_min, k_max
            k = vect_k[i]
            # if (bin1!=bin2):
            #     if (k-delta_k_int_1 > k_min ):
            #         extrem_1 = k- delta_k_int_1
            #     if (k+delta_k_int_1 < k_max ):
            #         extrem_2 = k + delta_k_int_1

            # Integration:
            params_int_k1[0] = k
            if i>5:
                average = np.average( integral[j_mu, i-5:i] )
            else:
                average = 42.
            if (k>0.07 and abs(average)<0.1):
                integral[j_mu,i]=0
            else:
                start = time.clock()
                integral[j_mu,i]  = sqrt(vol_shell_original(bin1)*vol_shell_original(bin2)) / (4*np.pi**2) * eval_integ_GSL_cquad(extrem_1, extrem_2, abs_prec_k1, rel_prec_k1, params_int_k1, W_k1, &F_k1, &n_eval)
                stop = time.clock()
                print "Bin: %d; (k,mu)=(%g, %g) --> %g --> %g sec." %(bin1,k,mu,integral[j_mu,i],stop-start)
            if (count%10==0):
                    np.savetxt(PATH_WINDOWED_SPECTRA+"%d_%d-var%d-mu%g_intAP.csv" %(bin1,bin2,num_var,mu), integral)
            count+=1
        np.savetxt(PATH_WINDOWED_SPECTRA+"%d_%d-var%d-mu%g_intAP.csv" %(bin1,bin2,num_var,mu), integral)
    np.savetxt(PATH_WINDOWED_SPECTRA+"%d_%d-var%d-mu%g_intAP.csv" %(bin1,bin2,num_var,mu), integral)
    return integral

def compute_AP_integral(var, bin):

    vect_k = np.logspace(np.log10(1e-3),np.log10(2e-1), 130)

    integral1 = AP_integral(var,bin,bin,vect_k)

def wrapper_intAP(var, bin1, bin2,kmin,kmax,Nk):
    vect_k=np.linspace(kmin,kmax,Nk)
    integral = AP_integral(var,bin1,bin2,vect_k)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(vect_k,integral,'r-',label="Int_DER")
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/int1/%d%d_var%d.pdf' %(bin1,bin2,var))
    return



# def plot_angDER(bin1,bin2,k,kp, N_z = 600,min_z=-1, max_z=1):
#     z_vect = np.linspace(min_z,max_z,N_z)
#     samples = np.zeros(N_z)
#     cdef double params[4]
#     params[0], params[1], params[2], params[3]= k, kp, bin1, bin2
#     for i in range(N_z):
#             samples[i] = arg_angle_int_derK(z_vect[i],params)

#     fig2=pl.figure()
#     ax1=fig2.add_subplot(111)
#     ax1.plot(z_vect,samples,'r-',label="Arg int. z --> DER K_ij")
#     # ax1.set_xlabel("$k$ [$h$/Mpc]")
#     # ax1.set_ylabel("[Mpc/$h$]^3")
#     ax1.grid(True)
#     ax1.legend(loc='best')
#     fig2.savefig('plots/intDER/arg_angle_intDER_%d%d_%g_%g.pdf' %(bin1,bin2,k,kp))
#     return

# def plot_intDER(bin1,bin2, k, N_kp = 600, kp_start=k_min, kp_stop=k_max, rel_input=rel_prec_mu1_DER, abs_input=abs_prec_mu1_DER):
#     global abs_prec_mu1_DER, rel_prec_mu1_DER
#     abs_prec_mu1_DER = abs_input
#     rel_prec_mu1_DER = rel_input
#     kp_vect = np.linspace(kp_start,kp_stop,N_kp)
#     samples = np.zeros(N_kp)
#     cdef double params[3]
#     params[0], params[1], params[2]= k, bin1, bin2
#     start = time.clock()
#     for i in range(N_kp):
#             samples[i] = arg_int_derK(kp_vect[i],params)
#     print time.clock()-start

#     fig2=pl.figure()
#     ax1=fig2.add_subplot(111)
#     ax1.plot(kp_vect,samples,'r-',label="Arg. int. kp --> DER K_ij")
#     ax1.grid(True)
#     ax1.legend(loc='best')
#     fig2.savefig('plots/intDER/arg_intDER_%d%d_%g.pdf' %(bin1,bin2,k))
#     return

# def plot_ang1(bin1,bin2,k,kp, num_var,  N_z = 600,min_z=-1, max_z=1):
#     z_vect = np.linspace(min_z,max_z,N_z)
#     samples = np.zeros(N_z)
#     cdef double params[4]
#     params[0], params[1], params[2], params[3]= k, kp, bin1, bin2
#     for i in range(N_z):
#             samples[i] = arg_integral_mu1(z_vect[i],params)

#     fig2=pl.figure()
#     ax1=fig2.add_subplot(111)
#     ax1.plot(z_vect,samples,'r-',label="Arg int. z --> Main term")
#     ax1.grid(True)
#     ax1.legend(loc='best')
#     fig2.savefig('plots/prove/arg_integral_mu1_%d%d_var%d_%g_%g.pdf' %(bin1,bin2,num_var,k,kp))
#     return


# def plot_int1(bin1,bin2, k, num_var, N_kp = 600, kp_start=k_min, kp_stop=k_max):
#     kp_vect = np.linspace(kp_start,kp_stop,N_kp)
#     samples = np.zeros(N_kp)
#     cdef double params[4]
#     params[0], params[1], params[2], params[3]= k, bin1, bin2, num_var
#     for i in range(N_kp):
#             samples[i] = arg_integral_k1(kp_vect[i],params)

#     fig2=pl.figure()
#     ax1=fig2.add_subplot(111)
#     ax1.plot(kp_vect,samples,'r-',label="Arg. int. kp --> Main term")
#     # ax1.set_xlabel("$k$ [$h$/Mpc]")
#     # ax1.set_ylabel("[Mpc/$h$]^3")
#     ax1.grid(True)
#     ax1.legend(loc='best')
#     fig2.savefig('plots/arg_integral_k1_%d%d_var%d_%g.pdf' %(bin1,bin2,num_var,k))
#     return

# # def plot_der_W1(bin1,bin2,k,kp, num_var, N_z = 600, min_z=-1, max_z=1):
# #     z_vect = np.linspace(min_z,max_z,N_z)
# #     samples = np.zeros(N_z)
# #     for i in range(N_z):
# #             samples[i] = der_W1(k,kp,bin1,bin2,z_vect[i])
# #     fig2=pl.figure()
# #     ax1=fig2.add_subplot(111)
# #     ax1.plot(z_vect,samples,'r-',label="der_W1")
# #     # ax1.set_xlabel("$k$ [$h$/Mpc]")
# #     # ax1.set_ylabel("[Mpc/$h$]^3")
# #     ax1.grid(True)
# #     ax1.legend(loc='best')
# #     fig2.savefig('plots/der_W1_%d%d_var%d_%g_%g.pdf' %(bin1,bin2,num_var,k,kp))
# #     return

# # def plot_der_Wcompiled(bin, N_kmod = 1000, start=0, end=0.5):
# #     k_vect = np.linspace(start,end,N_kmod)
# #     samples = np.zeros(N_kmod)
# #     for i in range(N_kmod):
# #             samples[i] = derW_compiled(k_vect[i],bin)

# #     fig2=pl.figure()
# #     ax1=fig2.add_subplot(111)
# #     ax1.plot(k_vect,samples,'r-',label="der_W_modK")
# #     # ax1.set_xlabel("$k$ [$h$/Mpc]")
# #     # ax1.set_ylabel("[Mpc/$h$]^3")
# #     ax1.grid(True)
# #     ax1.legend(loc='best')
# #     fig2.savefig('plots/der_W_modK.pdf' )
# #     return

# # def plot_Wcompiled(bin, N_kmod = 1000, start=0, end=0.5):
# #     k_vect = np.linspace(start,end,N_kmod)
# #     samples = np.zeros(N_kmod)
# #     for i in range(N_kmod):
# #             samples[i] = W_compiled(k_vect[i],bin)

# #     fig2=pl.figure()
# #     ax1=fig2.add_subplot(111)
# #     ax1.plot(k_vect,samples,'r-',label="der_W_modK")
# #     # ax1.set_xlabel("$k$ [$h$/Mpc]")
# #     # ax1.set_ylabel("[Mpc/$h$]^3")
# #     ax1.grid(True)
# #     ax1.legend(loc='best')
# #     fig2.savefig('plots/W_modK.pdf' )
# #     return

