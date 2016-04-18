import sys
import numpy as np
import sympy as sym
import scipy as scy
import matplotlib.pyplot as pl
import os
import time

INPUT_FOLDER = "OUTPUT/14bins_FFT/"
OUTPUT_FOLDER = "plots/marginalisation/comparison/14bins_FFT/"
name = "kmax0.2"
noAP = np.loadtxt(open(INPUT_FOLDER+"FM_uncorr_K0_"+name+".csv","rb"),delimiter=" ",skiprows=0)
yesAP = np.loadtxt(open(INPUT_FOLDER+"FM_uncorr_K1_"+name+".csv","rb"),delimiter=" ",skiprows=0)

# uFM_K1 = np.loadtxt(open("OUTPUT/FM_uncorr_K1.csv","rb"),delimiter=" ",skiprows=0)

# # Correction factor volume:
# correction_factor = np.pi * 15000. / (180.**2 * 4.)
# cFM, uFM = cFM*correction_factor, uFM*correction_factor

# Adding k_derivative term: (only non-correlated at the moment)
# uFM = uFM_K1

# ------------------------------------
# Checking dependence on k_max:
# ------------------------------------
# uFM_04 = np.loadtxt(open(INPUT_FOLDER+"../OUTPUT/FM_uncorr_0.45.csv","rb"),delimiter=" ",skiprows=0)
# uFM_02 = np.loadtxt(open(INPUT_FOLDER+"../OUTPUT/FM_uncorr_0.2.csv","rb"),delimiter=" ",skiprows=0)

# difference = cFM / uFM

# print abs((uFM_04-uFM_02) / uFM_02)

# np.savetxt(INPUT_FOLDER+"FM_difference.csv", difference )


# ------------------------------------
# Computing uncertaintes: (finally)
# ------------------------------------

# inv_cFM = np.linalg.inv(cFM)
# inv_uFM = np.linalg.inv(uFM)

n_var = {"h": 0, "n_s": 1, "\\Omega_b": 2, "\\Omega_c": 3, "w_0": 4, "w_1": 5, "\\sigma_8": 6, "b_0": 7,"b_1": 8,"b_2": 9,"b_3": 10,"b_4": 11,"b_5": 12,"b_6": 13, "b_7": 14}
var_names = ["h", "n_s", "\\Omega_b", "\\Omega_c", "w_0", "w_1", "\\sigma_8", "b_0","b_1","b_2","b_3","b_4","b_5","b_6", "b_7"]
ref_values = {'\\Omega_m': 0.3163288037424815,'\\Omega_c': 0.12/0.67**2,  'h': 0.67, '\\Omega_b':  0.022/0.67**2, 'n_s': 0.965, "w_0": -0.98, "w_1": 0., "\\sigma_8": 0.83, "b_0": 1.,"b_1": 1.,"b_2": 1.,"b_3": 1.,"b_4": 1.,"b_5": 1.,"b_6": 1., "b_7": 1.}

# Fully marginalised uncertaintes:
# for i in range(cFM.shape[0]):
#     sigmaC, sigmaU = np.sqrt(inv_cFM[i,i]), np.sqrt(inv_uFM[i,i])
#     print "%s\t\tFM_corr=%g \tFM_uncorr=%g  \t%g%%" %(var_names[i],sigmaC, sigmaU, (sigmaC-sigmaU)/sigmaC*100)


# Plot ellipse for two variables:

def marginalise(inv_matrix,var1,var2):
    temp_matrix = np.row_stack((inv_matrix[var1,:],inv_matrix[var2,:]))
    return np.column_stack((temp_matrix[:,var1], temp_matrix[:,var2]))

import cmath
# Return axes and angle:
def sigma1_ellipse(FM_matrix, var1, var2, scal=1.52):
    inv_matrix = np.linalg.inv(FM_matrix)
    sub_matrix = marginalise(inv_matrix,var1,var2)
    eigvals, eigvects = np.linalg.eig(sub_matrix)

    # Angle:
    for j in range(eigvects.shape[1]): # check columns
        if (eigvects[0,j]<0 and eigvects[1,j]<0) or (eigvects[0,j]>0 and eigvects[1,j]<0):
            eigvects[:,j] = -eigvects[:,j]
    u, v = eigvects[:,0], np.array([1,0])
    c = np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v) # -> cosine of the angle between
    angle = np.arccos(np.clip(c, -1, 1))

    # Axes:
    a, b = [scal*np.abs(cmath.sqrt(eigval)) for eigval in eigvals]

    # OLD METHOD:
    # # Find axes:
    # term1 = (sub_matrix[0,0]+sub_matrix[1,1])/2.
    # term2 = np.sqrt((sub_matrix[0,0]-sub_matrix[1,1])**2/4. +sub_matrix[0,1]**2)
    # a_2 = scal*np.sqrt( term1+term2 )
    # b_2 = scal*np.sqrt( abs(term1-term2) )
    # # eigenvalues = np.linalg.eig(inv_matrix)[0]

    # angle_2 = np.arctan( 2*sub_matrix[0,1]/(sub_matrix[0,0]-sub_matrix[1,1]) ) / 2.
    return 2*a, 2*b, angle



from math import degrees
import matplotlib.patches as patch

def print_ellipse(noAP,yesAP,var1,var2): #var names!
    if n_var[var1]>n_var[var2]:
        vartemp = var1
        var1 = var2
        var2 = vartemp
    results1 = sigma1_ellipse(noAP,n_var[var1],n_var[var2])
    results2 = sigma1_ellipse(yesAP,n_var[var1],n_var[var2])
    print results2[0], results2[1]
    print degrees(results1[2]), degrees(results2[2])
    ellipse1 = patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=results1[0], height=results1[1], angle=degrees(results1[2]), alpha= 0.3, color="g", label="No AP effect")
    # ellipse2 = patch.Ellipse(xy=(0,0), width=0.03, height=0.002, angle=46)
    ellipse2 = patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=results2[0], height=results2[1], angle=degrees(results2[2]), alpha= 0.3, color="b", label="With AP effect")
    fig=pl.figure()
    ax1=fig.add_subplot(111)
    ax1.add_patch(ellipse2)
    ax1.add_patch(ellipse1)

    # ax1.set_aspect('equal')
    if abs(results1[0]/2*np.cos(results1[2]))>abs(results1[1]/2*np.sin(results1[2])):
        x_range_1 = abs(results1[0]/2*np.cos(results1[2]))
    else:
        x_range_1 = abs(results1[1]/2*np.sin(results1[2]))
    if abs(results1[1]/2*np.cos(results1[2]))>abs(results1[0]/2*np.sin(results1[2])):
        y_range_1 = abs(results1[1]/2*np.cos(results1[2]))
    else:
        y_range_1 = abs(results1[0]/2*np.sin(results1[2]))
    if abs(results2[0]/2*np.cos(results2[2]))>abs(results2[1]/2*np.sin(results2[2])):
        x_range_2 = abs(results2[0]/2*np.cos(results2[2]))
    else:
        x_range_2 = abs(results2[1]/2*np.sin(results2[2]))
    if abs(results2[1]/2*np.cos(results2[2]))>abs(results2[0]/2*np.sin(results2[2])):
        y_range_2 = abs(results2[1]/2*np.cos(results2[2]))
    else:
        y_range_2 = abs(results2[0]/2*np.sin(results2[2]))

    x_range = x_range_1 if x_range_1>x_range_2 else x_range_2
    y_range = y_range_1 if y_range_1>y_range_2 else y_range_2


    # if results2[0]>results2[1]:
    #     x_range = results2[0]/2*np.cos(results2[2])
    #     y_range = results2[0]/2*np.sin(results2[2])
    # else:
    #     x_range = -results2[1]/2*np.sin(results2[2])
    #     y_range = results2[1]/2*np.cos(results2[2])
    # ax1.set_xlim(-1,1)
    # ax1.set_ylim(-1,1)
    ax1.set_xlim( ref_values[var1]-x_range*1.3,ref_values[var1]+x_range*1.3 )
    ax1.set_ylim( ref_values[var2]-y_range*1.3, ref_values[var2]+y_range*1.3 )
    ax1.set_xlabel("$%s$" %(var1))
    ax1.set_ylabel("$%s$" %(var2))
    ax1.grid(True)
    ax1.legend(loc='best')
    # ax1.legend().set_visible(False)
    fig.savefig(OUTPUT_FOLDER+'compareAP_'+name+'_vars_%d_%d.pdf' %(n_var[var1],n_var[var2]))



var1, var2 = "h", "\\Omega_b"
# var1, var2 = "h", "\\Omega_c"
# var1, var2 = "h", "n_s"
# var1, var2 = "n_s", "\\sigma_8"
# var1, var2 = "h", "w_0"
n_var1, n_var2 = n_var[var1], n_var[var2]
print_ellipse(noAP,yesAP,var1,var2)

# for i in range(15):
#     for j in range(i+1,15):
#         print_ellipse(uFM,cFM,var_names[i],var_names[j])
#
