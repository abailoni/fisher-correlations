
intAP_FOLDER = "INPUT/AP_integrals/2D_14bins/"

from scipy.interpolate import interp2d

AP_integral, AP_integral_corr = {}, {}
def import_AP_int():
    import_variables = ["Om_b","Om_c","w_0"]
    global AP_integral, AP_integral_corr
    #print "\nImporting and interpolating AP integrals..."
    for var in import_variables:
        AP_integral[var] = [0.]*N_bins
        AP_integral_corr[var] = [0.]*N_bins
        for bin in range(N_bins):
            #Please fix the integrals... and then take them out
            saved_data = vol_shell_original(bin)/gsl_pow_3(2*PI) / (sqrt(vol_shell_original(bin)*vol_shell_original(bin)) / (4*PI**2)) * np.loadtxt(open(intAP_FOLDER+"%d_%d-var%d-mu1_intAP.csv" %(bin,bin,n_var[var]),"rb"))
            vect_k = np.logspace(np.log10(1e-3),np.log10(2e-1), 130) #CHANGE!!!
            vect_mu = np.linspace(-1.,1.,5)
            AP_integral[var][bin] = interp2d(vect_k, vect_mu, saved_data, kind='cubic')
            if bin<N_bins-1:
                saved_data_corr = (4*np.pi**2)/gsl_pow_3(2*PI)  * np.loadtxt(open(intAP_FOLDER+"corr/%d_%d-var%d_intAP.csv" %(bin,bin+1,n_var[var]),"rb"))
                AP_integral_corr[var][bin] = interp2d(vect_k, vect_mu, saved_data_corr, kind='cubic')
    #print "--> Done!"


def AP_integral_fun(k,mu,bin,var):
    return AP_integral[var][bin](k,mu)
