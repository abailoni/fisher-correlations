import time
import sys, platform, os
from matplotlib import pyplot as plt
import numpy as np
#uncomment this if you are running remotely and want to keep in synch with repo changes
#if platform.system()!='Windows':
#    !cd $HOME/git/camb; git pull github master; git log -1
print('Using CAMB installed at '+ os.path.realpath(os.path.join(os.getcwd(),'..')))
sys.path.insert(0,os.path.realpath(os.path.join(os.getcwd(),'..')))
import camb
from camb import model, initialpower


h = 0.67
Omega_b = 0.049
Omega_c = 0.25

omegab_h2 = 0.022
omegac_h2 = 0.12
H0 = 67.

print Omega_b*h*h
print (omegab_h2)/(H0/100.)**2

c = 299792458.
print c*h/3000./1000.

#Now get matter power spectra and sigma8 at redshift 0 and 0.8
startTIME = time.time()
pars = camb.CAMBparams()
pars.set_cosmology(H0=h*100, ombh2=0.022, omch2=0.12)
# pars.set_cosmology()
pars.set_dark_energy(w=-0.98) #re-set defaults

pars.InitPower.set_params(ns=0.965)
#Not non-linear corrections couples to smaller scales than you want
pars.set_matter_power(redshifts=[0.], kmax=2.0)

#Linear spectra
pars.NonLinear = model.NonLinear_none

results = camb.get_results(pars)
kh, z, pk = results.get_matter_power_spectrum(minkh=0.02, maxkh=2, npoints = 200)
s8 = np.array(results.get_sigma8())
print z
print pk[0,0], s8
print pk[-1,0]*(0.83/s8[-1])**2

quit()
intermediate = time.time()
print intermediate-startTIME

# Change one paramter:
pars.set_cosmology(H0=70, ombh2=0.022, omch2=0.122)
results.calc_power_spectra(pars)
kh_2, z_2, pk_2 = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
s8_2 = np.array(results.get_sigma8())
print time.time()-intermediate


#Non-Linear spectra (Halofit)
pars.NonLinear = model.NonLinear_both
results.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

print results.get_sigma8()


for i, (redshift, line) in enumerate(zip(z,['-','--'])):
    plt.loglog(kh_2, pk_2[i,:], color='k', ls = line)
    plt.loglog(kh_nonlin, pk_nonlin[i,:], color='r', ls = line)
plt.xlabel('k/h Mpc');
plt.legend(['linear','non-linear'], loc='lower left');
# plt.title('Matter power at z=%s and z= %s'%tuple(z));
plt.show()
quit()
