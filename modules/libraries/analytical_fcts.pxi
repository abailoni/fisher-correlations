#**************************************************************
#**************************************************************
#
# Defining functions: (sympy --> internal cython)
#
#**************************************************************
#**************************************************************

#--------------------------------
# Start defining some variable:
#--------------------------------
Om_c, h, Om_b, w_0, w_1, n_s, gamma, sigma8 = sym.symbols('Om_c h Om_b w_0 w_1 n_s gamma sigma8')

# Other variables:
z = sym.symbols('z')
# Toy variables:
zx = sym.symbols('zx') #Actually this is always a number...

# Pivot redshift:
z_p = 0
# Cosmological constant w(z):
w = w_0 + w_1*(z-z_p)

# Dark Energy density (flat universe k=0):
Omega_DE = 1 - (Om_b+Om_c)


#-----------------------------------------
# Define Hubble parameter and distances:
#-----------------------------------------

# Hubble parameter: (diveded by H_0)
Hub = sym.sqrt( (Om_c+Om_b)*(1+z)**3 + Omega_DE*sym.exp(3* sym.integrate( (1+w.subs(z,zx))/(1+zx), (zx,0,z)) ) ) #[z, w_1, w_0, Om_m]

# Comoving distance in Mpc/h:
def comov_dist(zx,Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'],w_0=ref_values['w_0'],w_1=ref_values['w_1']):
    return c_H0*NInt(1/Hub,z,0,zx, w_1=w_1, Om_b=Om_b, Om_c=Om_c, w_0=w_0)

# Angular diameter distance: (in units of c/H_0)
def D_a(zx,Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'],w_0=ref_values['w_0'],w_1=ref_values['w_1']):
    return 1./(1+zx)*comov_dist(zx,Om_b,Om_c,w_0,w_1)


#--------------------------------------------------------------
# Growth factor and beta factor:
#--------------------------------------------------------------

Om_m_z = (Om_c+Om_b)* (1+z)**3 / Hub**2
#Om_m_z_fct = fnExpr(Om_m_z)
Om_m_z_py = SymToPy(Om_m_z)

def Growth(zx,Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'],gamma=ref_values['gamma'],w_1=ref_values['w_1'],w_0=ref_values['w_0']):
    return np.exp( NInt( Om_m_z**gamma/(1+z), z, zx, 0, w_1=w_1, Om_b=Om_b, Om_c=Om_c, gamma=gamma, w_0=w_0))

def beta(bin):
    return ( fnEv(Om_m_z_py,z=z_avg[bin],w_1=ref_values['w_1'],w_0=ref_values['w_0'],Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'])**ref_values['gamma'] / bias_bins[bin] )

#--------------------------------------------------------------
# Volume of a shell and top-hat function:
#--------------------------------------------------------------

# Partial volume of the sphere restriced to the survey area:
cdef double vol_shell(int bin):
    return PI*PI*survey_area/(180*180 * 3) * (-gsl_pow_3(com_zbin[bin]) + gsl_pow_3(com_zbin[bin+1]))
def vol_shell_py(bin):
    return vol_shell(bin)

# Normal 3D volume of a sphere:
cdef double vol_shell_original(int bin):
    return 4*PI/3. * (-gsl_pow_3(com_zbin[bin]) + gsl_pow_3(com_zbin[bin+1]))

def vol_shell_original_py(bin):
    return vol_shell_original(bin)

# Used only for an easier computation of W and K:
cdef double vol_shell_mod(int bin):
    return 1/3. * (-gsl_pow_3(com_zbin[bin]) + gsl_pow_3(com_zbin[bin+1]))

# Fourier transform of the top-hat function:
r, k, kp, r1, r2, r3, r4 = sym.symbols('r k kp r1 r2 r3 r4')
W_integral = sym.integrate( sym.sin(r*k)* r**2 / (k*r), (r,r1,r2)  ) #[k, r1, r2]
fun_W_integral = fnExpr(W_integral) #[k, r1, r2] --> bad optimized, not used

cdef double W_compiled(double mod_k, int bin):
    cdef double r1=com_zbin[bin], r2=com_zbin[bin+1]
    if mod_k!=0:
        return -(-r1*cos(mod_k*r1)/mod_k + sin(mod_k*r1)/(mod_k*mod_k))/mod_k + (-r2*cos(mod_k*r2)/mod_k + sin(mod_k*r2)/(mod_k*mod_k))/mod_k
    else:
        return 0.

#cdef double W(double k_modulus, int bin):
#    return 1./vol_shell_mod(bin)*W_compiled(k_modulus, com_zbin[bin],com_zbin[bin+1])

cdef double K(double mod_k, int bin1, int bin2):
    return 1./(vol_shell_mod(bin1)*vol_shell_mod(bin2)) * W_compiled(mod_k,bin1) * W_compiled(mod_k,bin2)
def K_py(mod_k,bin1,bin2):
    return K(mod_k,bin1,bin2)

# Compute k_max at each z: (actually here x_max is already k_max...)
def sigma_8(x_max, zx):
    integral_k = quad(lambda kp: kp**2 * zero_spectrum(kp) * Fourier_W_k(kp*np.pi/(2*x_max))**2, k_min, k_max, epsrel=1e-3)[0]
    integral_z = np.exp( -2. * NInt( Om_m_z**gamma/(1+z), z, 0., zx , gamma=ref_values['gamma'], w_1=ref_values['w_1'], w_0=ref_values['w_0'], Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'] ))
    return 1/(2.*PI**2) * integral_z * integral_k

def sigma_8_eq(x_max, zx):
    return np.sqrt(sigma_8(x_max, zx))-np.sqrt(sigma_8_equation)

def root_k_max(zx):
    return newton(sigma_8_eq,starting_point_sigma8,args=tuple([zx]))


#--------------------------------------------------------------
# Analytical derivatives for Om_m, w_0, w_1 and gamma:
#--------------------------------------------------------------

four_parameters = ['Om_b', 'Om_c', 'w_0', 'w_1','gamma']
mu, b_i, b_j = sym.symbols('mu b_i b_j')
redshift_factor = (1+Om_m_z**gamma/b_i*mu**2) * (1+Om_m_z**gamma/b_j*mu**2)

# Derivates of lnG and Beta wrt the four parameters:
# (for Beta, bias b_i added only to EUCLID data..)
lnG_der, Beta_der = {}, {}

for var in four_parameters:
    lnG_der[var] = lambda zx, Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'],gamma=ref_values['gamma'],w_1=ref_values['w_1'],w_0=ref_values['w_0'],var=var:  NInt(sym.diff(Om_m_z**sym.symbols('gamma'),sym.symbols(var))/(1+z), z, zx, 0., w_1=w_1, Om_b=Om_b, Om_c=Om_c, gamma=gamma, w_0=w_0)
    Beta_der[var] = lambda z, Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'],gamma=ref_values['gamma'],w_1=ref_values['w_1'],w_0=ref_values['w_0'], var=var: sym.diff(Om_m_z**sym.symbols('gamma'),sym.symbols(var)).subs([('w_1',w_1),('Om_b',Om_b),('Om_c',Om_c), ('gamma',gamma), ('w_0',w_0), ('z',z)])


# Derivatives of k and mu wrt ln(H) and ln(D):
cdef double mu_der_lnH(double mu):
    return(-mu * (mu**2-1))
cdef double mu_der_lnD(double mu):
    return(-mu * (mu**2-1))
cdef double k_der_lnH(double mu, double k):
    return(k * mu**2)
cdef double k_der_lnD(double mu, double k):
    return(k * (mu**2-1))

# Derivates of lnH and lnD wrt the four parameters:
lnH_der, lnD_der = {}, {}

for var in four_parameters: # num_var = [3-5] + gamma
    par = sym.symbols(var)
    lnH_der[var] = lambda z, Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'],w_1=ref_values['w_1'],w_0=ref_values['w_0'], par=par:  (sym.diff(Hub,par)/Hub).subs([('w_1',w_1),('Om_b',Om_b),('Om_c',Om_c), ('w_0',w_0), ('z',z)])
    lnD_der[var] = lambda z, Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'],w_1=ref_values['w_1'],w_0=ref_values['w_0'], var=var: 1./D_a(z,Om_b,Om_c,w_0,w_1) * 1./(1+z)*c_H0* (-1.) * quad(lambda zx: lnH_der[var](zx,Om_b,Om_c,w_1,w_0),0,z,epsrel=INT_PREC)[0]


# Derivative of mu wrt the four parameters: # num_var = [3-5] + gamma
cdef double mu_der(double mu, np.intp_t bin, np.intp_t var_num):
    return(mu_der_lnH(mu)*lnH_der_data[var_num][bin] + mu_der_lnD(mu)*lnD_der_data[var_num][bin])
cdef double k_der(double mu, double k, np.intp_t bin, np.intp_t var_num):
    return(k_der_lnH(mu,k)*lnH_der_data[var_num][bin] + k_der_lnD(mu,k)*lnD_der_data[var_num][bin])


#--------------------------------------------------------------
# Derivative of K in k:
#--------------------------------------------------------------

######################################
# IMPROVE WITH C pow()....!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
######################################

# Derivative wrt k of W:
cdef double derW_compiled(double mod_k, int bin):
    cdef double r1=com_zbin[bin], r2=com_zbin[bin+1]
    if mod_k!=0:
        return -(r1**2*sin(mod_k*r1)/mod_k + 2*r1*cos(mod_k*r1)/mod_k**2 - 2*sin(mod_k*r1)/mod_k**3)/mod_k + (r2**2*sin(mod_k*r2)/mod_k + 2*r2*cos(mod_k*r2)/mod_k**2 - 2*sin(mod_k*r2)/mod_k**3)/mod_k + (-r1*cos(mod_k*r1)/mod_k + sin(mod_k*r1)/mod_k**2)/mod_k**2 - (-r2*cos(mod_k*r2)/mod_k + sin(mod_k*r2)/mod_k**2)/mod_k**2
    else:
        return 0.


# Derivative of the square root wrt k:
cdef double sq_root_compiled(double k, double kp, double z):
    return (k - kp*z)/sqrt(k*k - 2*k*kp*z + kp*kp)


# W2(k) * der_W1(k):
cdef double W2_x_derW1(double k, double kp, int bin1, int bin2, double z):
    cdef:
        double sq_root = sqrt(k*k+kp*kp-2*k*kp*z)
        double der_sq_root
    if ( 1.-z< 1e-10 ): # Numerically z=1
        if ( kp>k ):
            der_sq_root = -1.
        else:
            der_sq_root = 1.
    else:
        if ( abs(k-kp)<5e-4 ): # kp is almost k
            der_sq_root = sqrt((1-z)/2.) #Simplified formula for the der.
            #print ".",
        else:
            der_sq_root = sq_root_compiled(k,kp,z)
    return 1./(vol_shell_mod(bin1)*vol_shell_mod(bin2)) * W_compiled(sq_root,bin2) * derW_compiled(sq_root,bin1) * der_sq_root

def wrapper_derW(k,bin1,bin2):
    return 1./(vol_shell_mod(bin1)*vol_shell_mod(bin2)) * W_compiled(k,bin2) * derW_compiled(k,bin1)

def W2_x_derW1_FFT_3D(vect_k,bin1,bin2):
    [kx, ky, kz] = np.meshgrid(vect_k,vect_k,vect_k)
    mod_k = np.sqrt(kx*kx + ky*ky + kz*kz)

    cdef int shape = vect_k.shape[0]
    results = np.empty((shape,shape,shape))
    cdef:
        double[::1] vect_k_c = vect_k
        double[:,:,::1] kz_c = kz
        double[:,:,::1] mod_k_c = mod_k
        double[:,:,::1] results_c = results
    # Can be optimised because it's symmetric:
    cdef double mod_k_val
    for i in range(shape):
        for j in range(shape):
            for t in range(shape):
                mod_k_val = mod_k_c[i,j,t]
                results_c[i,j,t]=1./(vol_shell_mod(bin1)*vol_shell_mod(bin2)) * W_compiled(mod_k_val,bin2) * derW_compiled(mod_k_val,bin1) * kz_c[i,j,t]/mod_k_val
    return results

# Some helpful external function for generating interpolated densities and
# bias from EUCLID data: (only z in [0.7, 2.0])
def bias_py(z):
    z_avg = np.array([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0])
    bias_bins = np.array([1.083, 1.125, 1.104, 1.126, 1.208, 1.243, 1.282, 1.292, 1.363, 1.497, 1.486, 1.491, 1.573, 1.568])
    return interp1d(z_avg,bias_bins,kind="slinear")(z)
def density_py(z): # EUCLID-2012
    z_avg = np.array([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0])
    n_dens = np.array([1.25, 1.92, 1.83, 1.68, 1.51, 1.35, 1.20, 1.00, 0.80, 0.58, 0.38, 0.35, 0.21, 0.11])
    return interp1d(z_avg,n_dens,kind="slinear")(z)


#**************************************************************
#**************************************************************
#
# Callable functions and storage variables:
#
#**************************************************************
#**************************************************************

#-----------------------------------------
# Cython memory views: (for storing data)
#-----------------------------------------
cdef:
    double[::1] com_zbin # N_bins+1
    double[::1] com_zbin_avg # N_bins
    #double[::1] DA_zbin = np.array([ D_a(zbn) for zbn in z_in])
    double[::1] Growth_bins # N_bins
    double[::1] beta_bins # N_bins
    double[:,::1] lnG_der_data # N_vars x N_bins
    double[:,::1] Beta_der_data # N_vars x N_bins
    double[:,::1] lnH_der_data # N_vars x N_bins
    double[:,::1] lnD_der_data # N_vars x N_bins
    double[::1] ref_val_v = np.empty(5) #for CLASS derivatives
    double[::1] k_max_data # N_bins
    double[::1] n_dens_c # N_bins

# For correlated FM:
cdef:
    double[:,::1] inverse_C_v
    double[:,::1] C_v
    double[:,::1] P_der_1_v
    double[:,::1] P_der_2_v
    double[:,::1] N_v
# Numpy matrices:
C = np.zeros([N_bins, N_bins])


#----------------------------------------------------
# Functions to compute distances and set new data:
#----------------------------------------------------

#*****************************
# set_ref_values(**args)
#*****************************
#
# Every time after updating the data it is necessary to compute
# again distancese and derivatives!
#
# Options:
#   - the names of the variables are: h, Om_m, Om_b,
#     n_s, gamma, Om_k, w_0, w_1

def set_ref_values(**args):
    var_names = ['h', 'Om_m', 'Om_b','Om_c', 'n_s', 'gamma', 'Om_k', 'w_0', 'w_1']
    for var_name in args:
        if var_name in var_names:
            ref_values[var_name] = args[var_name]
        else:
            print "Not valid option inserted!!"
            break

#*****************************
# set_survey(**options)
#*****************************
#
# Every time after updating the data it is necessary to compute
# again distancese and derivatives!
#
#
# Options:
#   - "bins_list": the python list should include both extremes.
#         Thus the last number is the border of the last bin
#   - "dens_list" is a python list containing only N_bins densities
#   - "survey_area": area in sq deg [default 15000]
#   - "bias_list" is a python list containting only N_bins biases

def set_survey(**args):
    options_val = ["bins_list", "survey_area", "dens_list", "bias_list"]
    for option in args:
        if option not in options_val:
            print "Not valid option inserted!!"
            break
        if option=="bins_list":
            global N_bins, z_in
            N_bins = len(args[option])-1
            z_in = np.array(args[option])
        if option=="dens_list":
            global n_dens
            n_dens = np.array(args[option])
        if option=="survey_area":
            global survey_area
            survey_area = args[option]
        if option=="bias_list":
            global bias_bins
            bias_bins = np.array(args[option])



#*****************************
# compute_survey_DATA()
#*****************************
# No input required
#
def compute_survey_DATA():
    # Setting bins:
    global z_avg, dz
    z_avg = np.array([ (z_in[i]+z_in[i+1])/2. for i in range(N_bins)])
    dz = np.array([z_in[i+1]-z_in[i] for i in range(N_bins)])
    # Distances:
    print "\nComputing survey data:"
    print " - distances..."
    global com_zbin, com_zbin_avg
    com_zbin = np.array([ comov_dist(zbn) for zbn in z_in])
    com_zbin_avg = np.array([ (com_zbin[i]+com_zbin[i+1])/2. for i in range(N_bins)])
    # Derivatives and funtions:
    print " - derivatives..."
    global lnG_der_data, Beta_der_data, lnH_der_data, lnD_der_data, Growth_bins, beta_bins
    Growth_bins = np.array([Growth(zx) for zx in z_avg])
    beta_bins = np.array([ beta(bin) for bin in range(N_bins)])
    lnG_der_data, Beta_der_data, lnH_der_data, lnD_der_data = np.zeros([N_vars,N_bins]), np.zeros([N_vars,N_bins]), np.zeros([N_vars,N_bins]), np.zeros([N_vars,N_bins])
    #four_parameters = ['Om_m', 'w_0', 'w_1', 'gamma']
    four_parameters = ['Om_b','Om_c', 'w_0', 'w_1']
    for var in four_parameters: # num_var = [3-5] + gamma
        for bin in range(N_bins):
            # for the first two also add bias:
            lnG_der_data[n_var[var]][bin] = lnG_der[var](z_avg[bin])
            Beta_der_data[n_var[var]][bin] = 1./bias_bins[bin] * Beta_der[var](z_avg[bin])
            lnH_der_data[n_var[var]][bin] = lnH_der[var](z_avg[bin])
            lnD_der_data[n_var[var]][bin] = lnD_der[var](z_avg[bin])
    ## Test and print data:
    #for bin in range(N_bins):
    #    print lnH_der_data[3][bin], lnH_der_data[4][bin], lnH_der_data[5][bin], lnH_der_data[3][bin]/lnH_der_data[5][bin]
    # Other data: (for CLASS derivatives)
    ref_values_arr = [0, ref_values['h'], ref_values['n_s'], ref_values['Om_b'], ref_values['Om_c']]
    for i in range(5):
        ref_val_v[i] = epsilon * ref_values_arr[i]
    # Densities:
    global n_dens_c
    n_dens_c = efficiency*n_dens*1e-3

    ## k_max:
    #print " - k_max at each z..."
    #global k_max_data
    #k_max_data = np.array([root_k_max(zx) for zx in z_avg])

    # Correlation FM data:
    global N_tot_vars
    N_tot_vars = 7 + N_bins
    global P_der_1_v, P_der_2_v, N_v, C, C_v
    P_der_1_v, P_der_2_v = np.zeros([N_bins, N_bins]), np.zeros([N_bins, N_bins])
    inv_n_dens = 1./(n_dens*efficiency*1e-3)
    N_v = np.identity(N_bins) * inv_n_dens
    C = np.zeros([N_bins, N_bins])
    C_v = C
    print "--> Done!\n"
    return
