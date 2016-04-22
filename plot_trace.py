import sys
import numpy as np
import sympy as sym
import scipy as scy
import matplotlib
import matplotlib.pyplot as pl
import os
import time


from scipy.interpolate import interp1d
import modules.cFM as cFM

cFM.init()
k_vect, trace = cFM.plot_trace(0,0,50,10,1e-3,0.2)
prova = interp1d(k_vect,trace,kind="cubic")
k_vect_2 = np.linspace(1e-3,0.2,3000)
trace_int = prova(k_vect_2)

fig2=pl.figure()
ax1=fig2.add_subplot(111)
ax1.plot(k_vect,trace,'r-',label="Trace along $\\mu=-1$")
ax1.plot(k_vect_2,trace_int,'b-',label="Trace interp $\\mu=-1$")
# ax1.plot(vect_k,class_fct['P_0'](vect_k),'b-',label="Class Spectrum")
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_yscale('log')
# ax1.set_xlabel("$k$ [$h$/Mpc]")
# ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
fig2.savefig('plots/trace/vabbe.pdf')
# np.savetxt("OUTPUT/trace/trace_along_mu_vars%d%d.csv"%(var1,var2),np.column_stack((k_vect,samples[:,0])))


# cFM.plot_trace(1,1,100,10,1e-3,0.2)
