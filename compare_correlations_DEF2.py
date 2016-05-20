import sys
import numpy as np
import sympy as sym
import scipy as scy
import matplotlib
import matplotlib.pyplot as pl
import os
import time

# n_var = {"h": 0, "n_s": 1, "\\Omega_b": 2, "\\Omega_c": 3, "w_0": 4, "w_1": 5, "\\sigma_8": 6, "b_0": 7,"b_1": 8,"b_2": 9,"b_3": 10,"b_4": 11,"b_5": 12,"b_6": 13, "b_7": 14}
# var_names = ["h", "n_s", "\\Omega_b", "\\Omega_c", "w_0", "w_1", "\\sigma_8"]
n_var = {"h": 0, "n_s": 1, "\\Omega_b": 2, "\\Omega_c": 3, "w_0": 4, "\\sigma_8": 5, "b_0": 6,"b_1": 7,"b_2": 8,"b_3": 9,"b_4": 10,"b_5": 11,"b_6": 12, "b_7": 13}
var_names = ["h", "n_s", "\\Omega_b", "\\Omega_c", "w_0", "\\sigma_8"]

ref_values = {'\\Omega_m': 0.3163288037424815,'\\Omega_c': 0.12/0.67**2,  'h': 0.67, '\\Omega_b':  0.022/0.67**2, 'n_s': 0.965, "w_0": -0.98, "w_1": 0., "\\sigma_8": 0.83, "b_0": 1.,"b_1": 1.,"b_2": 1.,"b_3": 1.,"b_4": 1.,"b_5": 1.,"b_6": 1., "b_7": 1.}


class FishMatr(object):
    def __init__(self):
        return 0.



# -----------------------------------------------
# ROUNTINES:
# -------------------------------------------


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
    # return 2*a, 2*b, angle
    # Invert axes:
    return 2*a, 2*b, np.pi/2. - angle



from math import degrees
import matplotlib.patches as patch


# Computes the ellipses wrt var1 and var2 for all the FM passed.
#
# Inputs:
#   - matrices is a list of matrices
# The output is:
#   - a list of ellipses
#   - a list of [ x_range, y_range ]

def compute_ellipses(matrices, var1, var2): #var names
    if n_var[var1]>n_var[var2]:
        vartemp = var1
        var1 = var2
        var2 = vartemp
    N = len(matrices)
    ellipses, ranges = [], np.empty((N,2))
    for matrix, xyRange in zip(matrices, ranges):
        results = sigma1_ellipse(matrix,n_var[var1],n_var[var2])
        ellipses.append(patch.Ellipse(xy=(ref_values[var2], ref_values[var1]), width=results[0], height=results[1], angle=degrees(results[2])))
        if abs(results[0]/2*np.cos(results[2]))>abs(results[1]/2*np.sin(results[2])):
            xyRange[0] = abs(results[0]/2*np.cos(results[2]))
        else:
            xyRange[0] = abs(results[1]/2*np.sin(results[2]))
        if abs(results[1]/2*np.cos(results[2]))>abs(results[0]/2*np.sin(results[2])):
            xyRange[1] = abs(results[1]/2*np.cos(results[2]))
        else:
            xyRange[1] = abs(results[0]/2*np.sin(results[2]))
    return ellipses, ranges


def max_range(ranges):
    ranges = np.array(ranges)
    if len(ranges.shape)==2:
        return np.amax(ranges, axis=0)
    elif len(ranges.shape)==3: # all vars
        return np.amax(ranges, axis=1)


def print_ellipses(ellipses, maxRange, var1, var2, pdf_name, labels, legend="best",colors=['r','b','g'], opacities=[0.3]*10):
    if n_var[var1]<n_var[var2]: # exchange variables axes
        vartemp = var1
        var1 = var2
        var2 = vartemp

    N = len(ellipses)
    fig=pl.figure()
    ax1=fig.add_subplot(111)
    for ellipse, label, color, opacity in zip(ellipses, labels, colors, opacities):
        ellipse.set_label(label)
        ellipse.set_alpha(opacity)
        ellipse.set_color(color)
        ax1.add_patch(ellipse)
    ax1.set_xlim( ref_values[var1]-maxRange[0]*1.3,ref_values[var1]+maxRange[0]*1.3 )
    ax1.set_ylim( ref_values[var2]-maxRange[1]*1.3, ref_values[var2]+maxRange[1]*1.3 )
    ax1.set_xlabel("$%s$" %(var1))
    ax1.set_ylabel("$%s$" %(var2))
    ax1.grid(True)
    ax1.legend(loc=legend)
    # ax1.legend().set_visible(False)
    fig.savefig(OUTPUT_FOLDER+pdf_name+'.pdf')


def compute_ellipses_allVars(matrices, N_vars=len(var_names)): #var names
    N = len(matrices)
    all_ellipses, all_ranges = [], []
    for nvar2 in range(N_vars):
        for nvar1 in range(nvar2+1,N_vars):
            var1, var2 = var_names[nvar1], var_names[nvar2]
            ellipses, ranges = [], np.empty((N,2))
            for matrix, xyRange in zip(matrices, ranges):
                results = sigma1_ellipse(matrix,n_var[var1],n_var[var2])
                ellipses.append(patch.Ellipse(xy=(ref_values[var2], ref_values[var1]), width=results[0], height=results[1], angle=degrees(results[2])))
                if abs(results[0]/2*np.cos(results[2]))>abs(results[1]/2*np.sin(results[2])):
                    xyRange[0] = abs(results[0]/2*np.cos(results[2]))
                else:
                    xyRange[0] = abs(results[1]/2*np.sin(results[2]))
                if abs(results[1]/2*np.cos(results[2]))>abs(results[0]/2*np.sin(results[2])):
                    xyRange[1] = abs(results[1]/2*np.cos(results[2]))
                else:
                    xyRange[1] = abs(results[0]/2*np.sin(results[2]))
            all_ellipses = all_ellipses + ellipses
            all_ranges = all_ranges + [ranges]
    return all_ellipses, all_ranges


def print_ellipses_allVars(all_ellipses, maxRanges, pdf_name, labels, legend="best",colors=['r','b','g'], N_vars=len(var_names),ticks_size=3, opacities=[0.3]*10):
    N_plots = N_vars*N_vars - N_vars*(N_vars+1)/2
    N_ellipses = len(all_ellipses)
    # FM per plot:
    N_FM = N_ellipses/N_plots

    font_ticks = {'weight' : 'normal',
        'size'   : ticks_size}
    matplotlib.rc('font', **font_ticks)

    # Reshape the ellipses:
    all_ellipses = zip(*[iter(all_ellipses)]*N_FM)
    fig=pl.figure()
    ell_index = 0
    for nvar2 in range(N_vars):
        for nvar1 in range(nvar2+1,N_vars):
            var2, var1 = var_names[nvar1], var_names[nvar2] #exchange axes
            plot_index = (nvar1-1)*(N_vars-1) + (nvar2+1)
            ax=fig.add_subplot(N_vars-1,N_vars-1,plot_index)
            for ellipse, label, color, opacity in zip( list(all_ellipses[ell_index]), labels, colors, opacities):
                ellipse.set_label(label)
                ellipse.set_alpha(opacity)
                ellipse.set_color(color)
                ax.add_patch(ellipse)
            ax.set_xlim( ref_values[var1]-maxRange[ell_index,0]*1.3,ref_values[var1]+maxRange[ell_index,0]*1.3 )
            ax.set_ylim( ref_values[var2]-maxRange[ell_index,1]*1.3, ref_values[var2]+maxRange[ell_index,1]*1.3 )
            if nvar1==N_vars-1:
                ax.set_xlabel("$%s$" %(var1), fontsize=15)
            if nvar2==0:
                ax.set_ylabel("$%s$" %(var2), fontsize=15)
            ax.grid(True)
            ax.legend(loc=legend)
            ax.legend().set_visible(False)
            ell_index+=1
    fig.set_tight_layout(True)
    # fig.tight_layout(pad=0.1, w_pad=0.1, h_pad=0.1)
    fig.savefig(OUTPUT_FOLDER+pdf_name+'.pdf')

# -------------------------------------------
# END
# -------------------------------------------


INPUT_FOLDER = "OUTPUT/DEF_AP0/"
OUTPUT_FOLDER = "plots/marginalisation/comparison/DEF_AP0/"
name = "EUCLID2012_DEF_AP0-14bins"
uFM = np.loadtxt(open(INPUT_FOLDER+"FMcorr_"+name+"-uncorr.csv","rb"),delimiter=" ",skiprows=0)
wFM = np.loadtxt(open(INPUT_FOLDER+"FMcorr_"+name+"-windowFun.csv","rb"),delimiter=" ",skiprows=0)
cFM = np.loadtxt(open(INPUT_FOLDER+"FMcorr_"+name+"-correlations+windowFun.csv","rb"),delimiter=" ",skiprows=0)




var1, var2 = "h", "\\Omega_b"
# var1, var2 = "\\Omega_b", "\\Omega_c"
# var1, var2 = "h", "n_s"
# var1, var2 = "n_s", "\\sigma_8"
# var1, var2 = "h", "w_0"
n_var1, n_var2 = n_var[var1], n_var[var2]

ellipses, ranges = compute_ellipses_allVars([cFM,wFM,uFM],3)
maxRange = max_range(ranges)
print_ellipses_allVars(ellipses,maxRange,"EUCLID2012_AP0-14bins",["Top hat + correlations","Top hat","Standard"],"upper left",['g','b','r'],3,9)

# ellipses, ranges = compute_ellipses([wFM,cFM,uFM],var1,var2)
# maxRange = max_range(ranges)
# print_ellipses(ellipses,maxRange,var1,var2,"compare",["Top hat","Top hat + correlations","Standard"],"upper left",['g','b','r'])

