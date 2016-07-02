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
Hub_data = SymToLambda(Hub,numpy=True,**ref_values)


# Comoving distance in Mpc/h:
#def comov_dist(zx,Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'],w_0=ref_values['w_0'],w_1=ref_values['w_1']):
#    return c_H0*NInt(1/Hub,z,0,zx, w_1=w_1, Om_b=Om_b, Om_c=Om_c, w_0=w_0)

def Hubble(z,**cosmo_params):
    return Lambda_Ev(Hub_data,z,**cosmo_params)

# Comoving distance:
def comov_dist(z,**cosmo_params):
    # If there are no options use the default values:
    for param in ref_values:
        if param not in cosmo_params:
            cosmo_params[param] = ref_values[param]
    # Check if z is an array:
    starting_point = 0. if not numpy_check(z) else np.zeros(z.shape)
    return c_H0 * NIntegrate(1/Hub,'z',starting_point,z,1e-6,**cosmo_params)

# Angular diameter distance: (in units of c/H_0)
def D_a(z,**cosmo_params):
    return 1./(1+z)*comov_dist(z,**cosmo_params)


#--------------------------------------------------------------
# Growth factor and beta factor:
#--------------------------------------------------------------

Om_m_z = (Om_c+Om_b)* (1+z)**3 / Hub**2
#Om_m_z_fct = fnExpr(Om_m_z)
#Om_m_z_py = SymToPy(Om_m_z)
Om_m_z_data = SymToLambda(Om_m_z,numpy=True,**ref_values)


def Growth(zx,**cosmo_params):
    for param in ref_values:
        if param not in cosmo_params:
            cosmo_params[param] = ref_values[param]
    # Check if zx is an array:
    starting_point = 0. if not numpy_check(zx) else np.zeros(zx.shape)
    return np.exp( NIntegrate(Om_m_z**gamma/(1+z), 'z', zx, starting_point, 1e-6, **cosmo_params))

def beta(bins,**cosmo_params):
    """ bins can be a vector """
    if 'gamma' not in cosmo_params:
        cosmo_params['gamma']=ref_values['gamma']
    return  Lambda_Ev(Om_m_z_data,z_avg[bins],**cosmo_params)**cosmo_params['gamma'] / bias_bins_numpy[bins]

# Just for test:
def growth_rate_f(bin,**cosmo_params):
    """ The input 'bin' can be a vector """
    if 'gamma' not in cosmo_params:
        cosmo_params['gamma']=ref_values['gamma']
    return Lambda_Ev(Om_m_z_data,z_avg[bin],**cosmo_params)**cosmo_params['gamma']

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

def vol_shell_original_py(bins):
    """ The input 'bins' can be a vector """
    if numpy_check(bins):
        results = np.empty(bins.shape)
        for i, bin in enumerate(bins):
            results[i] = vol_shell_original(bin)
        return results
    else:
        return vol_shell_original(bins)

# Used for the numerical derivatives:
def vol_spherical_shell(bins,**cosmo_params):
    """ The input 'bins' can be a vector """
    if not numpy_check(bins):
        bins = np.array([bins])
    bins_indices = np.append(bins,[bins.max()+1])
    com_dists = comov_dist(z_in[bins_indices],**cosmo_params)
    return 4*PI/3. * (-np.power(com_dists[:-1],3) + np.power(com_dists[1:],3))

# Used only for an easier computation of W and K:
cdef double vol_shell_mod(int bin):
    return 1/3. * (-gsl_pow_3(com_zbin[bin]) + gsl_pow_3(com_zbin[bin+1]))

# Fourier transform of the top-hat function:
r, k, kp, r1, r2, r3, r4 = sym.symbols('r k kp r1 r2 r3 r4')
W_integral = sym.integrate( sym.sin(r*k)* r**2 / (k*r), (r,r1,r2),conds='none')

# Lambdify optimized to be used in a vectorized way. See K_FFTconvolution()
W_integral_data = SymToLambda(W_integral,numpy=True,order_vars=['k'],r1=1.,r2=1.)

cdef double W_compiled(double mod_k, int bin):
    cdef double r1=com_zbin[bin], r2=com_zbin[bin+1]
    if mod_k!=0:
        return -(-r1*cos(mod_k*r1)/mod_k + sin(mod_k*r1)/(mod_k*mod_k))/mod_k + (-r2*cos(mod_k*r2)/mod_k + sin(mod_k*r2)/(mod_k*mod_k))/mod_k
    else:
        print "Warning: W_fourier(k) computed at k=0"
        return vol_shell_mod(bin)

def W_py(mod_k,bin):
    return W_compiled(mod_k,bin)

#def W_py_2(mod_k,bin):
#    return Lambda_Ev(W_integral_data, mod_k, r1=com_zbin[bin], r2=com_zbin[bin+1])

# bin_pairs is a matrix of dim (2, N_pairs)
# One column of 'allPairs_comDist' contains (R_i, R_j, R_i+1, R_j+1). Every colum represent a bin_pair (i,j). Same for 'allPairs_comVol' but there just (i,j).
# All the constants are already included. CHANGE OTHER FUNCTIONS!!!!! --> integral_1
# Results is of the shape (N_pairs, N_k)

def K_FFTconvolution(bin_pairs, mods_k, **cosmo_params):

    # Compute comov. distances and volumes:
    bins = np.arange(bin_pairs.min(),bin_pairs.max()+1)
    com_distances = comov_dist(z_in[np.concatenate((bins,[bins[-1]+1]))],**cosmo_params)
    com_volumes = vol_spherical_shell(bins,**cosmo_params)
    # For each pair:
    allPairs_bins = np.row_stack((bin_pairs, bin_pairs+1))
    allPairs_comDist, allPairs_comVol = com_distances[allPairs_bins], com_volumes[bin_pairs]

    # Check for zero mods_k:
    idxs_k, idxs_bins = (mods_k!=0).nonzero(), np.arange(bin_pairs.shape[1])
    idxMat_k, idxMat_bin = np.meshgrid(idxs_k,idxs_bins)

    # Compute K(mods_k) for each pair:
    allPairs_mods_k = np.tile(mods_k, (bin_pairs.shape[1],1))

    #print "W:"
    #print Lambda_Ev(W_integral_data, mods_k, r1=com_distances[0], r2=com_distances[1])[:10]
    #print Lambda_Ev(W_integral_data, allPairs_mods_k[idxMat_bin,idxMat_k], r1=allPairs_comDist[0], r2=allPairs_comDist[2])[0,:10]
    #print "done"
    result_matrix = Lambda_Ev(W_integral_data, allPairs_mods_k[idxMat_bin,idxMat_k], r1=allPairs_comDist[0], r2=allPairs_comDist[2]) * Lambda_Ev(W_integral_data, allPairs_mods_k[idxMat_bin,idxMat_k], r1=allPairs_comDist[1], r2=allPairs_comDist[3])
    results = (2./ (PI*np.sqrt(np.prod(allPairs_comVol,axis=0))) * result_matrix.T).T

    # Adjust zero mods_k:
    idxs_zerok = tuple((mods_k==0).nonzero()[0])
    results = np.insert(results.T, idxs_zerok, np.sqrt(np.prod(allPairs_comVol,axis=0))/gsl_pow_3(2*PI), axis=0).T

    return results




#cdef double W(double k_modulus, int bin):
#    return 1./vol_shell_mod(bin)*W_compiled(k_modulus, com_zbin[bin],com_zbin[bin+1])

cdef double K(double mod_k, int bin1, int bin2):
    if mod_k>1e-6:
        return 1./(vol_shell_mod(bin1)*vol_shell_mod(bin2)) * W_compiled(mod_k,bin1) * W_compiled(mod_k,bin2)
    else:
        return 1.
def K_py(mod_k,bin1,bin2):
    return K(mod_k,bin1,bin2)

# Compute k_max at each z: (actually here x_max is already k_max...)
def sigma_8(x_max, zx, **cosmo_params):
    for param in ref_values:
        if param not in cosmo_params:
            cosmo_params[param] = ref_values[param]
    integral_k = quad(lambda kp: kp**2 * zero_spectrum(kp) * Fourier_W_k(kp*np.pi/(2*x_max))**2, k_min, k_max, epsrel=1e-3)[0]
    integral_z = np.exp( -2. * NInt( Om_m_z**gamma/(1+z),z,0.,zx,**cosmo_params))
    return 1/(2.*PI**2) * integral_z * integral_k

def sigma_8_eq(x_max, zx, **cosmo_params):
    return np.sqrt(sigma_8(x_max, zx, **cosmo_params))-np.sqrt(sigma_8_equation)

# MISSING cosmo_params.....!!!!!!!!
def root_k_max(zx):
    return newton(sigma_8_eq,starting_point_sigma8,args=tuple([zx]))


#--------------------------------------------------------------
# Analytical derivatives for Om_m, w_0, w_1 and gamma:
#--------------------------------------------------------------

#parameters_derivated = ['Om_b', 'Om_c', 'w_0', 'w_1','gamma']
parameters_derivated = ['Om_b', 'Om_c', 'w_0']
mu, b_i, b_j = sym.symbols('mu b_i b_j')
redshift_factor = (1+Om_m_z**gamma/b_i*mu**2) * (1+Om_m_z**gamma/b_j*mu**2)

# Derivates of lnG and Beta wrt the four parameters:
# (for Beta, bias b_i added only to EUCLID data.. Why? This is actually not Beta but the growth rate f)

#lnG_der, Beta_der = {}, {}
#for var in parameters_derivated:
#    lnG_der[var] = lambda zx,var=var:  NInt(sym.diff(Om_m_z**sym.symbols('gamma'),sym.symbols(var))/(1+z), z, zx, 0., **ref_values)
#    Beta_der[var] = lambda z,var=var: sym.diff(Om_m_z**sym.symbols('gamma'),sym.symbols(var)).subs([('w_1',ref_values['w_1']),('Om_b',ref_values['Om_b']),('Om_c',ref_values['Om_c']), ('gamma',ref_values['gamma']), ('w_0',ref_values['w_0']), ('z',z)])

beta_der_Lambda_data = {}
for var in parameters_derivated:
    beta_der_Lambda_data[var] = SymToLambda(sym.diff(Om_m_z**sym.symbols('gamma'),sym.symbols(var)),numpy=True,**ref_values)

lnG_der_lamda = {}
for var in parameters_derivated:
    lnG_der_lamda[var] = integrationFun(sym.diff(Om_m_z**sym.symbols('gamma'),sym.symbols(var))/(1+z),'z',**ref_values)

# Analytical ones: (alwasy computed in default values...)
def lnG_der(zx,var):
    # Check if zx is an array:
    starting_point = 0. if not numpy_check(zx) else np.zeros(zx.shape)
    return NIntegrate_fun(lnG_der_lamda[var], zx, starting_point,1e-6)
def Beta_der(z,var):
    return Lambda_Ev(beta_der_Lambda_data[var],z)


# Derivatives of k and mu wrt ln(H) and ln(D):
cdef double mu_der_lnH(double mu):
    return(-mu * (mu**2-1))
cdef double mu_der_lnD(double mu):
    return(-mu * (mu**2-1))
cdef double k_der_lnH(double mu, double k):
    return(k * mu**2)
cdef double k_der_lnD(double mu, double k):
    return(k * (mu**2-1))

#lnH_der, lnD_der = {}, {}
#for var in parameters_derivated: # num_var = [3-5] + gamma
#    par = sym.symbols(var)
#    lnH_der[var] = lambda z, Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'],w_1=ref_values['w_1'],w_0=ref_values['w_0'], par=par:  (sym.diff(Hub,par)/Hub).subs([('w_1',w_1),('Om_b',Om_b),('Om_c',Om_c), ('w_0',w_0), ('z',z)])
#    lnD_der[var] = lambda z, Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c'],w_1=ref_values['w_1'],w_0=ref_values['w_0'], var=var: 1./D_a(z,Om_b,Om_c,w_0,w_1) * 1./(1+z)*c_H0* (-1.) * quad(lambda zx: lnH_der[var](zx,Om_b,Om_c,w_1,w_0)/fnEv(Hub_Theano,z=zx,w_1=ref_values['w_1'],w_0=ref_values['w_0'],Om_b=ref_values['Om_b'],Om_c=ref_values['Om_c']) ,0,z,epsrel=INT_PREC)[0]

lnH_der_Lambda_data = {}
for var in parameters_derivated:
    lnH_der_Lambda_data[var] = SymToLambda(sym.diff(Hub,sym.symbols(var))/Hub,numpy=True,**ref_values)


# Derivates of lnH and lnD wrt the four parameters: (analytical)
def lnH_der(z,var):
    return Lambda_Ev(lnH_der_Lambda_data[var],z)
# This is a mess...!!!
def lnD_der(z,var):
    return 1./D_a(z) * 1./(1+z)*c_H0* (-1.) * quad(lambda zx: lnH_der(zx,var)/fnEv(Hub_py,z=zx,**ref_values),0,z,epsrel=INT_PREC)[0]


# NUMERICAL DERIVATIVES:
def Fun_der_num(fun,z,var): # for G, H and D_a
    """ The input z can be a vector """
    return (fun(z,**{var: ref_values[var]+epsilon}) - fun(z,**{var: ref_values[var]-epsilon}) ) / (2*epsilon)
def Fun_der_num_bins(fun,bin, var): # for beta
    """ The input bin can be a vector """
    return (fun(bin,**{var: ref_values[var]+epsilon}) - fun(bin,**{var: ref_values[var]-epsilon}) ) / (2*epsilon)
def Fun_der_redshift(fun,z,**cosmo_params):  # for G, H and D_a
    """ The input z can be a vector """
    return (fun(z+epsilon,**cosmo_params) - fun(z-epsilon,**cosmo_params)) / (2*epsilon)



# Derivative of mu wrt the four parameters: # num_var = [3-5] + gamma
cdef double mu_der(double mu, np.intp_t bin, np.intp_t var_num):
    return mu_der_lnH(mu)*lnH_der_data[var_num][bin] + mu_der_lnD(mu)*lnD_der_data[var_num][bin]
cdef double k_der(double mu, double k, np.intp_t bin, np.intp_t var_num):
    return k_der_lnH(mu,k)*lnH_der_data[var_num][bin] + k_der_lnD(mu,k)*lnD_der_data[var_num][bin]

def k_der_py(mu,k,bin,var_num):
    return k_der(mu,k,bin,var_num)



# Checking AP numer. derivative:
cdef double AP_fact_R(double mu, int bin, int nvar): # n_var [2-4]
    return sqrt( (Hub_mod[nvar+1][bin]*Da_mod[nvar+1][bin]*mu)**2 - (Hub_mod[0][bin]*Da_mod[0][bin])**2 *(mu*mu-1) ) / (Hub_mod[0][bin]*Da_mod[nvar+1][bin])
cdef double k_AP(double k, double mu, int bin, int nvar): # n_var [2-4]
    return k * AP_fact_R(mu,bin,nvar)
cdef double mu_AP(double k, double mu, int bin, int nvar): # n_var [2-4]
    return mu * Hub_mod[nvar+1][bin] / (AP_fact_R(mu,bin,nvar)*Hub_mod[0][bin])
cdef double Pk_AP_num_der(double k, double mu, int bin, int nvar): #n_var [2-4]
    cdef double kAP = k_AP(k,mu,bin,nvar)
    return ( zero_spectrum(kAP)-zero_spectrum(k) )/ref_val_v[nvar+1]
# Mu AP term:
cdef double mu_num_term(double k, double mu, int bin, int nvar): #n_var [2-4]
    cdef double muAP = mu_AP(k,mu,bin,nvar)
    return 2 * (log(1+beta_mod[nvar+1][bin]*muAP*muAP) - log(1+beta_bins[bin]*mu*mu)) / ref_val_v[nvar+1]

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
# INVENTED DATA AT z=0.65 and z=2.05:
def bias_py(z):
    z_avg = np.array([0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.05])
    bias_bins = np.array([1.05, 1.083, 1.125, 1.104, 1.126, 1.208, 1.243, 1.282, 1.292, 1.363, 1.497, 1.486, 1.491, 1.573, 1.568, 1.58,])
    return interp1d(z_avg,bias_bins,kind="slinear")(z)
def density_py(z): # EUCLID-2012
    z_avg = np.array([0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5])
    n_dens = np.array([0.95 ,1.25, 1.92, 1.83, 1.68, 1.51, 1.35, 1.20, 1.00, 0.80, 0.58, 0.38, 0.35, 0.21, 0.11, 0.05])
    return interp1d(z_avg,n_dens,kind="slinear")(z)

# IMPROVED VERSION: it doesn't just interpolate at z_avg, but with the z_in vector as input it does an average in the bin:
def bias_advanced_py(z_in):
    # Interpolation:
    z_avg_default = np.array([0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5])
    bias_bins_default = np.array([1.05, 1.083, 1.125, 1.104, 1.126, 1.208, 1.243, 1.282, 1.292, 1.363, 1.497, 1.486, 1.491, 1.573, 1.568, 1.58,])
    interpolation = interp1d(z_avg_default,bias_bins_default,kind="slinear")

    # Average values:
    N_samples = 100
    z_samples = np.linspace(0.65,2.05,N_samples)
    bias_samples = interpolation(z_samples)
    results = np.empty(z_in.shape[0]-1)
    for i, z_stop in enumerate(z_in[1:]):
        indxs = np.logical_and(z_samples>=z_in[i], z_samples<=z_stop).nonzero()
        results[i] = np.average(bias_samples[indxs])
    return results


def density_advanced_py(z_in): # EUCLID-2012
    # Interpolation:
    z_avg_default = np.array([0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5])
    n_dens_default = np.array([0.95 ,1.25, 1.92, 1.83, 1.68, 1.51, 1.35, 1.20, 1.00, 0.80, 0.58, 0.38, 0.35, 0.21, 0.11, 0.05])
    interpolation = interp1d(z_avg_default,n_dens_default,kind="slinear")

    # Average values:
    N_samples = 100
    z_samples = np.linspace(0.65,2.05,N_samples)
    nDens_samples = interpolation(z_samples)
    results = np.empty(z_in.shape[0]-1)
    for i, z_stop in enumerate(z_in[1:]):
        indxs = np.logical_and(z_samples>=z_in[i], z_samples<=z_stop).nonzero()
        results[i] = np.average(nDens_samples[indxs])
    return results


# AP integral:
# Sympy variables:
phi1Sy, k1Sy, mu1Sy, muSy, lnH_Sy, lnD_Sy  = sym.symbols(' phi1Sy k1Sy mu1Sy muSy lnH_Sy lnD_Sy')

old_mu_squared = (1-muSy*muSy)*(1-mu1Sy*mu1Sy)*sym.cos(phi1Sy)**2 + (muSy*mu1Sy)**2 + 2 * (muSy*sym.sqrt(1-muSy**2)) * (mu1Sy*sym.sqrt(1-mu1Sy**2)) * sym.cos(phi1Sy)


k_der_sym = (lnH_Sy * k1Sy * old_mu_squared ) + (lnD_Sy * k1Sy * (old_mu_squared-1.))

phi1_intregral = sym.integrate(k_der_sym, (phi1Sy, 0., 2*PI))
phi1_intregral_data = SymToLambda(phi1_intregral,order_vars=["k1Sy", "mu1Sy", "muSy", "lnH_Sy", "lnD_Sy"])



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
    double[:,::1] lnVolShells_der_data # N_vars x N_bins
    #double[::1] lnH_der_0 # N_vars
    #double[::1] lnD_der_0 # N_vars
    double[::1] ref_val_v = np.empty(20) #for CLASS derivatives
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
sqrt_volume_shells = np.zeros([N_bins, N_bins])
P_der_1, P_der_2 = np.zeros([N_bins, N_bins]), np.zeros([N_bins, N_bins])

# Check AP numer. derivative:
cdef:
    double[:,::1] Hub_mod # N_vars x N_bins
    double[:,::1] Da_mod # N_vars x N_bins
    double[:,::1] beta_mod # N_vars x N_bins


#----------------------------------------------------
# Functions to compute distances and set new data:
#----------------------------------------------------

#*****************************
# set_FM_vars(**args)
#*****************************
#

def set_FM_vars(vars=available_FM_vars):
    global FM_vars_numbers, N_cosm_vars
    N_cosm_vars = len(vars)
    FM_vars_numbers = [0]*len(vars)
    for i, var in enumerate(vars):
        FM_vars_numbers[i] = n_FM_vars[var]

set_FM_vars() #This should be put in the init function...

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
            global bias_bins, bias_bins_numpy
            bias_bins_numpy = np.array(args[option])
            bias_bins = bias_bins_numpy


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
    com_zbin = comov_dist(z_in)
    com_zbin_avg = np.array([ (com_zbin[i]+com_zbin[i+1])/2. for i in range(N_bins)])

    # Derivatives and funtions:
    print " - derivatives..."
    global lnG_der_data, Beta_der_data, lnH_der_data, lnD_der_data, Growth_bins, beta_bins, lnVolShells_der_data
    Growth_numpy = Growth(z_avg)
    Growth_bins = Growth_numpy # Cython memoryview
    beta_bins = beta(np.arange(N_bins))
    lnG_der_numpy, Beta_der_numpy, lnH_der_numpy, lnD_der_numpy, lnVolShells_der_numpy = np.zeros([N_vars,N_bins]), np.zeros([N_vars,N_bins]), np.zeros([N_vars,N_bins]), np.zeros([N_vars,N_bins]), np.zeros([N_vars,N_bins])

    parameters_derivated = ['Om_b','Om_c', 'w_0']
    for var in parameters_derivated: # num_var = [2-4]
        # NUMERICAL:
        #Beta_der_numpy[n_var[var],:] = Fun_der_num_bins(beta,range(N_bins),var)
        #lnG_der_numpy[n_var[var],:] = Fun_der_num(Growth,z_avg,var) / Growth_numpy
        #lnH_der_numpy[n_var[var],:] = Fun_der_num(Hubble,z_avg,var) /Hubble(z_avg)
        lnD_der_numpy[n_var[var],:] = Fun_der_num(D_a,z_avg,var) / D_a(z_avg)
        lnVolShells_der_numpy[n_var[var],:] = Fun_der_num_bins(vol_spherical_shell,np.arange(N_bins),var) / vol_shell_original_py(np.arange(N_bins))

        # ANALYTICAL:
        lnG_der_numpy[n_var[var],:] = lnG_der(z_avg,var)
        Beta_der_numpy[n_var[var],:] = 1./bias_bins_numpy * Beta_der(z_avg,var)
        lnH_der_numpy[n_var[var],:] = lnH_der(z_avg,var)
        #lnD_der_numpy[n_var[var],:] = lnD_der(z_avg,var)

    # Save faster cython memoryvies:
    lnG_der_data, Beta_der_data, lnH_der_data, lnD_der_data, lnVolShells_der_data = lnG_der_numpy, Beta_der_numpy, lnH_der_numpy, lnD_der_numpy, lnVolShells_der_numpy

    # Other data: (for CLASS derivatives)
    ref_values_arr = [0, ref_values['h'], ref_values['n_s'], ref_values['Om_b'], ref_values['Om_c'], ref_values['w_0'], ref_values['w_1']]
    for i, val in enumerate(ref_values_arr):
        ref_val_v[i] = epsilon * val

    # Densities:
    global n_dens_c
    n_dens_c = efficiency*n_dens*1e-3

    ## k_max:
    #print " - k_max at each z..."
    #global k_max_data
    #k_max_data = np.array([root_k_max(zx) for zx in z_avg])

    # Correlation matrices initialisation:
    print " - correlation's matrices initialisation..."
    global sqrt_volume_shells
    #N_tot_vars = N_cosm_vars + N_bins
    sqrt_volume_shells = np.zeros([N_bins, N_bins])
    for bin1 in range(N_bins):
        for bin2 in range(bin1,N_bins):
            sqrt_volume_shells[bin1,bin2] = sqrt(sqrt(vol_shell(bin1)*vol_shell(bin2)))
            sqrt_volume_shells[bin2,bin1] = sqrt_volume_shells[bin1,bin2]

    global P_der_1, P_der_1_v, P_der_2, P_der_2_v, N_v, C, C_v
    P_der_1, P_der_2 = np.zeros([N_bins, N_bins]), np.zeros([N_bins, N_bins])
    P_der_1_v, P_der_2_v = P_der_1, P_der_2
    inv_n_dens = 1./(n_dens*efficiency*1e-3)
    N_v = np.identity(N_bins) * inv_n_dens
    C = np.zeros([N_bins, N_bins])
    C_v = C

    ## Checking AP numer. derivative:
    #print " - Test AP..."
    #global Hub_mod, Da_mod, beta_mod
    #Hub_mod, Da_mod, beta_mod = np.zeros([N_cosm_vars, N_bins]), np.zeros([N_cosm_vars, N_bins]), np.zeros([N_cosm_vars, N_bins])
    #dict_vars = {'Om_b': ref_values['Om_b'], 'Om_c': ref_values['Om_c'], 'w_0': ref_values['w_0'], 'w_1': ref_values['w_1']}
    #indices_again = {'Om_b':0, 'Om_c':1, 'w_0':2}
    #arr_vars = [ref_values['Om_b'], ref_values['Om_c'], ref_values['w_0']]
    #parameters_der = ['spectrum', 'Om_b','Om_c', 'w_0']
    #for var in parameters_der:
    #    for bin in range(N_bins):
    #        if var=="spectrum":
    #            Hub_mod[0][bin] = fnEv(Hub_Theano,z=z_avg[bin],**dict_vars)
    #            Da_mod[0][bin]  = D_a(z_avg[bin],*arr_vars)
    #        else:
    #            nvar = n_var_import[var]
    #            dict_temp, arr_temp = dict_vars.copy(), list(arr_vars)
    #            dict_temp[var], arr_temp[indices_again[var]] = dict_temp[var]+ref_val_v[nvar], arr_temp[indices_again[var]]+ref_val_v[nvar]
    #            Hub_mod[nvar][bin] = fnEv(Hub_Theano,z=z_avg[bin],**dict_temp)
    #            Da_mod[nvar][bin]  = D_a(z_avg[bin],*arr_temp)
    #            beta_mod[nvar][bin]= beta(bin,*arr_temp)
    return

# Export der_variables for Santi-comparison:
def print_DER():
    var_names = ['Om_b', 'Om_c', 'w_0']
    Hub = np.array([[lnH_der_data[n_var[var]][bin] for bin in range(N_bins)] for var in var_names])
    G = np.array([[lnG_der_data[n_var[var]][bin] for bin in range(N_bins)] for var in var_names])
    Da = np.array([[lnD_der_data[n_var[var]][bin] for bin in range(N_bins)] for var in var_names])
    Beta = np.array([[Beta_der_data[n_var[var]][bin] for bin in range(N_bins)] for var in var_names])
    return [Hub, G, Da, Beta]

