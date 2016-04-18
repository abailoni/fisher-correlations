import sys
import numpy as np
import sympy as sym
import scipy as scy
import matplotlib.pyplot as pl
import os
import time

from scipy.integrate import quad
from scipy.optimize import newton
from os import walk
from re import match
from sympy import srepr
from re import findall


import modules.uFM as uFM

uFM.compute_CAMB_spectra()

N_kp = 1000
kp_vect = np.linspace(0.0001,0.99,N_kp)
samples = np.zeros(N_kp)
for i in range(N_kp):
    samples[i] = uFM.CLASS_der_py(kp_vect[i],1)

fig2=pl.figure()
ax1=fig2.add_subplot(111)
ax1.plot(kp_vect,samples,'r-',label="Arg. int. kp --> Main term")
# ax1.set_xlabel("$k$ [$h$/Mpc]")
# ax1.set_ylabel("[Mpc/$h$]^3")
# ax1.set_xscale("log")
# ax1.set_yscale("log")
ax1.grid(True)
ax1.legend(loc='best')
fig2.savefig('plots/test_CAMB_spectrum.pdf')
