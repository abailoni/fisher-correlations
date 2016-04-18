# ###########################
# PLOT OBSERVED SPECTRUM:
# ###########################
def plot_observed_spectrum(bin1,bin2,N_k=500, N_mu=20):
    # N_k, N_mu = int(sys.argv[1]), int(sys.argv[2])
    # time_vect = np.zeros([N_k*N_mu])
    # time_inverse = np.zeros([N_k*N_mu])
    k_vect = np.linspace(k_min+0.001,k_max-0.04,N_k)
    mu_vect = np.linspace(-1.,1.,N_mu)
    samples = np.zeros([N_k,N_mu])
    for i_k in range(N_k):
        for i_mu in range(N_mu):
            samples[i_k,i_mu] = observed_spectrum(bin1,bin2,k_vect[i_k],mu_vect[i_mu])

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter

    fig = pl.figure()
    ax = fig.gca(projection='3d')
    mu_vect_m, k_vect_m  = np.meshgrid(mu_vect, k_vect)
    surf = ax.plot_surface(k_vect_m, mu_vect_m, samples, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    # ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    fig.savefig('plots/observed_spectrum_3D_bins%d%d.pdf'%(bin1,bin2))
    np.savetxt("output/observed_spectrum_bins%d%d.csv" %(bin1,bin2),samples)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(k_vect,samples[:,0],'r-',label="Observed spectrum along $\\mu=-1$")
    ax1.grid(True)
    ax1.legend(loc='best')
    # ax1.set_yscale('log')
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
    fig2.savefig('plots/obs_spctr_along_mu_bins%d%d.pdf'%(bin1,bin2))
    np.savetxt("output/obs_spctr_along_mu_bins%d%d.csv"%(bin1,bin2),np.column_stack((k_vect,samples[:,0])))
    return



# ###########################
# PLOT derivative_type_B:
# ###########################
def plot_der_B(bin1,bin2,var,double k_m=k_min+0.001, double k_M=k_max-0.04, N_k=1000, N_mu=10):
    k_vect = np.linspace(k_m,k_M,N_k)
    mu_vect = np.linspace(-1.,1.,N_mu)
    samples = np.zeros([N_k,N_mu])
    for i_k in range(N_k):
        for i_mu in range(N_mu):
            samples[i_k,i_mu] = der_type_B(bin1, bin2, k_vect[i_k], mu_vect[i_mu], var)

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter

    fig = pl.figure()
    ax = fig.gca(projection='3d')
    mu_vect_m, k_vect_m  = np.meshgrid(mu_vect, k_vect)
    surf = ax.plot_surface(k_vect_m, mu_vect_m, samples, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    # ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    fig.savefig('plots/der_typeB_3D_bins%d%d_var_%d.pdf'%(bin1,bin2,var))
    np.savetxt("output/der_typeB_bins%d%d_var_%d.csv" %(bin1,bin2,var),samples)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(k_vect,samples[:,0],'r-',label="Der. type B for $\\mu=-1$")
    ax1.set_xlabel("$k$ [$h$/Mpc]")
    ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/der_typeB_along_mu_bins%d%d_var%d.pdf'%(bin1,bin2,var))
    np.savetxt("output/der_typeB_along_mu_bins%d%d_var%d.csv"%(bin1,bin2,var),np.column_stack((k_vect,samples[:,0])))
    return

# ###########################
# PLOT derivative_type_A:
# ###########################
def plot_der_A(bin1,bin2,var,double k_m=k_min+0.001, double k_M=k_max-0.04, N_k=1000, N_mu=10):
    k_vect = np.linspace(k_m,k_M,N_k)
    mu_vect = np.linspace(-1.,1.,N_mu)
    samples = np.zeros([N_k,N_mu])
    for i_k in range(N_k):
        for i_mu in range(N_mu):
            samples[i_k,i_mu] = der_type_A(bin1,bin2,k_vect[i_k],-1.,var)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(k_vect,samples[:,0],'r-',label="Der. type A for $\\mu=-1$")
    ax1.set_xlabel("$k$ [$h$/Mpc]")
    ax1.set_ylabel("[Mpc/$h$]^3")
    ax1.grid(True)
    ax1.legend(loc='best')
    fig2.savefig('plots/der_typeA_along_mu_bins%d%d_var%d.pdf'%(bin1,bin2,var))
    np.savetxt("output/der_typeA_along_mu_bins%d%d_var%d.csv"%(bin1,bin2,var),np.column_stack((k_vect,samples[:,0])))
    return

# #################
# PLOTTING TRACE:
# #################

def plot_trace(var1,var2,N_k=1000,N_mu=10):
    k_vect = np.linspace(k_min+0.001,k_max-0.04,N_k)
    mu_vect = np.linspace(-1.,1.,N_mu)
    samples = np.zeros([N_k,N_mu])
    time_count = 0
    total_start = time.clock()
    for i_k in range(N_k):
        for i_mu in range(N_mu):
            samples[i_k,i_mu] = trace(k_vect[i_k],mu_vect[i_mu],var1,var2)
    total_stop = time.clock()
    np.savetxt("output/trace_3D_vars_%d%d.csv" %(var1,var2),samples)

    ######################
    # Plot this thing...
    # Result ---> it's a mess ;D
    ######################
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    fig = pl.figure()
    ax = fig.gca(projection='3d')
    # X = np.arange(-5, 5, 0.25)
    # Y = np.arange(-5, 5, 0.25)
    # R = np.sqrt(X**2 + Y**2)
    # Z = np.sin(R)
    mu_vect_m, k_vect_m  = np.meshgrid(mu_vect, k_vect)
    surf = ax.plot_surface(k_vect_m, mu_vect_m, samples, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    # # ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    fig.savefig('plots/trace_3D_vars_%d%d.pdf'%(var1,var2))

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(k_vect,samples[:,0],'r-',label="Trace along $\\mu=-1$")
    # ax1.plot(vect_k,class_fct['P_0'](vect_k),'b-',label="Class Spectrum")
    ax1.grid(True)
    ax1.legend(loc='best')
    # ax1.set_yscale('log')
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
    fig2.savefig('plots/trace_along_mu_vars%d%d.pdf'%(var1,var2))
    np.savetxt("output/trace_along_mu_vars%d%d.csv"%(var1,var2),np.column_stack((k_vect,samples[:,0])))
    return



# ###########################
# Plot argument k integral:
# ###########################
def FM_does_not_work(int var1, int var2, double k_min=new_k_min, double k_max=new_k_max, int N_k=1000):
    # N_k, N_mu = 1000, 10
    k_vect = np.linspace(k_min,k_max,N_k)
    samples = np.empty(N_k)
    # F_k.function = &argument_k
    # F_mu.function = &argument_mu
    cdef double params[2]
    # sdfsdf[0] = 3.
    params[0], params[1] = var1, var2
    for i_k in range(N_k):
        samples[i_k] = argument_k(k_vect[i_k],params)

    fig2=pl.figure()
    ax1=fig2.add_subplot(111)
    ax1.plot(k_vect,samples,'r-',label="Argument $k$ integral ($\\mu$ integrated)")
    ax1.grid(True)
    ax1.legend(loc='best')
    # ax1.set_yscale('log')
    # ax1.set_xlabel("$k$ [$h$/Mpc]")
    # ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
    fig2.savefig('plots/FM_arg_k_integral_vars%d%d.pdf' %(var1,var2))
    np.savetxt("output/FM_arg_k_integral_vars%d%d.csv" %(var1,var2),np.column_stack((k_vect,samples)))
    return
