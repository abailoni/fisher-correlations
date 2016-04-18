

# Importing modules and defining default data:

include "libraries/header.pxi"

include "libraries/import_CLASS.pxi"

include "libraries/CAMB.pxi"

include "libraries/analytical_fcts.pxi"


# Not necessary:
def all():
    compute_CAMB_spectra()
    compute_survey_DATA()

#---------------------------------------
# COMPUTING INTEGRAL 1:
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
            P_samples[i] = CLASS_der(vect_k[i],n_var_import[name_var])

    # Radial FFT:
    return np.sqrt(vol_shell_original(bin1)*vol_shell_original(bin2))/(2*np.pi)**3 *FFTt.radial_convolution(P_samples,K_samples,R)


# Storing and interpolating int1 with GSL:
cdef enum:
    max_N_bins = 50

cdef:
    interpolation_tools integral_1_tools[max_N_bins][max_N_bins][N_vars]
    # interpolation_tools integral_DER[max_N_bins][max_N_bins]

def store_int1():
    print "Computing, storing and interpolating windowed data (int1)..."
    vect_k = np.linspace(0.0001,0.5,5000)
    import_variables = ["spectrum","h","Om_b","Om_m","n_s"]
    start = time.time()
    for bin1 in range(N_bins):
        for bin2 in range(bin1,N_bins):
            print " - (%d,%d) --> " %(bin1,bin2),
            for name_var in import_variables:
                alloc_interp_GSL(vect_k, integral_1(bin1,bin2,name_var,vect_k), &integral_1_tools[bin1][bin2][n_var_import[name_var]])
    print "Done! (%g sec.)" %(time.time()-start)

