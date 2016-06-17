import sys
import numpy as np
import sympy as sym
import scipy as scy
import matplotlib
import matplotlib.pyplot as pl
import os
import time
import matplotlib.gridspec as gridspec


# n_var = {"h": 0, "n_s": 1, "\\Omega_b": 2, "\\Omega_c": 3, "w_0": 4, "w_1": 5, "\\sigma_8": 6, "b_0": 7,"b_1": 8,"b_2": 9,"b_3": 10,"b_4": 11,"b_5": 12,"b_6": 13, "b_7": 14}
# var_names = ["h", "n_s", "\\Omega_b", "\\Omega_c", "w_0", "w_1", "\\sigma_8"]
n_var = {"h": 0, "n_s": 1, "\\Omega^{(0)}_{\\mathrm{b}}": 2, "\\Omega^{(0)}_{\\mathrm{cdm}}": 3, "w_0": 4, "\\sigma_8": 5, "b_0": 6,"b_1": 7,"b_2": 8,"b_3": 9,"b_4": 10,"b_5": 11,"b_6": 12, "b_7": 13}
var_names = ["h", "n_s", "\\Omega^{(0)}_{\\mathrm{b}}", "\\Omega^{(0)}_{\\mathrm{cdm}}", "w_0", "\\sigma_8"]

ref_values = {'\\Omega_m': 0.3163288037424815,'\\Omega^{(0)}_{\\mathrm{cdm}}': 0.12/0.67**2,  'h': 0.67, '\\Omega^{(0)}_{\\mathrm{b}}':  0.022/0.67**2, 'n_s': 0.965, "w_0": -0.98, "w_1": 0., "\\sigma_8": 0.83, "b_0": 1.,"b_1": 1.,"b_2": 1.,"b_3": 1.,"b_4": 1.,"b_5": 1.,"b_6": 1., "b_7": 1.}


class FishMatr(object):
    def __init__(self):
        return 0.



# -----------------------------------------------
# ROUNTINES:
# -------------------------------------------


# BE CAREFUL: here it returns the inverse of the marginalised FM...
# Put variables to keep! (array)
def marginalise(FM_matrix,num_vars):
    inv_matrix = np.linalg.inv(FM_matrix)
    temp_matrix = np.take(inv_matrix,indices=num_vars,axis=1)
    return np.take(temp_matrix,indices=num_vars,axis=0)

# Put variables to maximise! (array)
def maximise(FM_matrix,num_vars_delete):
    temp_matrix = np.delete(FM_matrix,num_vars_delete,axis=1)
    return np.delete(temp_matrix,num_vars_delete,axis=0)

# Full marginalized uncertainties:
def uncertainty(FM_matrix,var):
    return np.sqrt(marginalise(FM_matrix,n_var[var]))

import cmath
# Return axes and angle:
def sigma1_ellipse(FM_matrix, num_var1, num_var2, scal=1.52):
    sub_matrix = marginalise(FM_matrix,[num_var1,num_var2])
    if num_var1>num_var2:
        sub_matrix = np.rot90(sub_matrix, k=2)

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

    # # OLD METHOD:
    # # Find axes:
    # term1 = (sub_matrix[0,0]+sub_matrix[1,1])/2.
    # term2 = np.sqrt((sub_matrix[0,0]-sub_matrix[1,1])**2/4. +sub_matrix[0,1]**2)
    # a_2 = scal*np.sqrt( term1+term2 )
    # b_2 = scal*np.sqrt( term1-term2 )
    # # eigenvalues = np.linalg.eig(inv_matrix)[0]

    # angle_2 = np.arctan( 2*sub_matrix[0,1]/(sub_matrix[0,0]-sub_matrix[1,1]) ) / 2.
    # return 2*a_2, 2*b_2, angle_2

    if num_var1>num_var2:
        return 2*a, 2*b, np.pi/2. - angle
    else:
        return 2*b, 2*a, -(np.pi/2. - angle)



from math import degrees
import matplotlib.patches as patch


def check_matplot_arguments(type_of_plot,**kargs):
    if not "labels" in kargs:
        kargs["labels"] = [""]*20
    if not "legend" in kargs:
        kargs["legend"] = "best"
    if not "colors" in kargs:
        kargs["colors"] = ['r','b','g','y','m','k']*20

    # For ellipses:
    if type_of_plot=="ellipses":
        if not "ticks_size" in kargs:
            kargs["ticks_size"] = 9
        if not "label_size" in kargs:
            kargs["label_size"] = kargs["ticks_size"]
        if not "opacities" in kargs:
            kargs["opacities"] = [0.3]*10

    # For line plots (particularly spectra...)
    elif type_of_plot=="linePlot":
        if not "ticks_size" in kargs:
            kargs["ticks_size"] = 13
        if not "label_size" in kargs:
            kargs["label_size"] = kargs["ticks_size"]
        if not "opacities" in kargs:
            kargs["opacities"] = [1.]*20
        if not "log" in kargs:
            kargs["log"] = ""
        if not "lineStyles" in kargs:
            kargs["lineStyles"] = ['-']*20
        if not "xyLabels" in kargs:
            kargs["xyLabels"] = ['$k$ [$h$/Mpc]', '$P(k)$ [(Mpc/$h$)$^3$]']
        if not "xrange" in kargs:
            kargs["xrange"] = 0
        if not "yrange" in kargs:
            kargs["yrange"] = 0
        if not 'grid' in kargs:
            kargs['grid'] = True
    return kargs

# -----------------------------------
# ELLIPSES:
# -----------------------------------

# Return the matplot ellipses wrt var1 and var2 for all the FM passed.
#
# Inputs:
#   - matrices is a list of matrices
# The output is:
#   - a list of ellipses
#   - [x_range, y_range] (in order to include all the considered ellipses)
def matplot_ellipses(matrices, var1, var2): #var names
    N = len(matrices)
    ellipses, ranges = [], np.empty((N,2))
    for i, matrix in enumerate(matrices):
        results = sigma1_ellipse(matrix,n_var[var1],n_var[var2])
        ellipses.append(patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=results[0], height=results[1], angle=degrees(results[2])))
        if abs(results[0]/2*np.cos(results[2]))>abs(results[1]/2*np.sin(results[2])):
            ranges[i,0] = abs(results[0]/2*np.cos(results[2]))
        else:
            ranges[i,0] = abs(results[1]/2*np.sin(results[2]))
        if abs(results[1]/2*np.cos(results[2]))>abs(results[0]/2*np.sin(results[2])):
            ranges[i,1] = abs(results[1]/2*np.cos(results[2]))
        else:
            ranges[i,1] = abs(results[0]/2*np.sin(results[2]))
    return ellipses, np.amax(ranges, axis=0)


# def matplot_ellipses_moreVars(matrices, vars=var_names): #var names
#     all_ellipses, all_ranges, all_varCouples = [], [], []
#     for var1 in vars[:-1]:
#         for var2 in vars[1:]:
#             ellipses, ranges = matplot_ellipses(matrices, var1, var2)
#             all_ellipses.append(ellipses), all_ranges.append(ranges)
#             all_varCouples.append([var1,var2])
#     return all_ellipses, all_ranges, all_varCouples


# def max_range(ranges):
#     ranges = np.array(ranges)
#     if len(ranges.shape)==2:
#         return np.amax(ranges, axis=0)
#     elif len(ranges.shape)==3: # all vars
#         return np.amax(ranges, axis=1)


# Plot ellipses on an axis:
def ellipses_onePlot(axis, matrices, var1, var2, **plot_kargs):
    plot_kargs = check_matplot_arguments("ellipses",**plot_kargs)

    matplotlib.rcParams['text.usetex'] = True
    font_style = {'weight' : 'normal', 'size': plot_kargs['ticks_size'],'family':'serif','serif':['Palatino']}
    matplotlib.rc('font', **font_style)

    ellipses_data = matplot_ellipses(matrices, var1, var2)
    for ellipse, label, color, opacity in zip(ellipses_data[0], plot_kargs['labels'], plot_kargs['colors'], plot_kargs['opacities']):
        ellipse.set_label(label)
        ellipse.set_alpha(opacity)
        ellipse.set_color(color)
        axis.add_patch(ellipse)

    ranges = ellipses_data[1]
    axis.set_xlim( ref_values[var1]-ranges[0]*1.3,ref_values[var1]+ranges[0]*1.3 )
    axis.set_ylim( ref_values[var2]-ranges[1]*1.3, ref_values[var2]+ranges[1]*1.3 )
    if "xyLabels" not in plot_kargs:
        axis.set_xlabel("$%s$" %(var1), size=plot_kargs["label_size"])
        axis.set_ylabel("$%s$" %(var2), size=plot_kargs["label_size"])
    else:
        axis.set_xlabel("$%s$" %(plot_kargs["xyLabels"][0]), size=plot_kargs["label_size"])
        axis.set_ylabel("$%s$" %(plot_kargs["xyLabels"][1]), size=plot_kargs["label_size"])
    axis.grid(True)
    axis.legend(loc=plot_kargs['legend'])

    return axis


# It works also with only two variables
# The x offset is before, the y one is after...
def ellipses_varsGrid(matrices, vars=var_names, x_offsetGrid=0, y_offsetGrid=0, **plot_kargs):
    Nvars = len(vars)

    # Create grid:
    gs = gridspec.GridSpec(x_offsetGrid+Nvars-1, y_offsetGrid+Nvars-1)
    axes = []

    # Print ellipses:
    plot_kargs["xyLabels"] = ["",""]
    for j, var1 in enumerate(vars[:-1]):
        for i, var2 in zip(range(j,Nvars),vars[1+j:]):
            ax_ij = pl.subplot(gs[x_offsetGrid+i, j])
            plot_kargs['xyLabels'][0] = var1 if i==Nvars-2 else ""
            plot_kargs['xyLabels'][1] = var2 if j==0 else ""
            axes.append( ellipses_onePlot(ax_ij, matrices, var1, var2, **plot_kargs) )
    return gs



# -----------------------------------
# 1D GAUSSIAN:
# -----------------------------------

def gaussian(x, mu, sigma):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sigma, 2.)))

def gaussian_matplot(axis, matrices, var, **plot_kargs):
    # Compute uncertainty:
    sigmas = np.array([ uncertainty(matrix,var) for matrix in matrices])
    mu = ref_values[var]
    x_vect = np.linspace(mu-np.max(sigmas)*2, mu+np.max(sigmas)*2,3000)
    ys = [gaussian(x_vect,mu,sigma) for sigma in sigmas]

    # Plot:
    plot_kargs["yrange"] = [0,1.1]
    if "xyLabels" not in plot_kargs:
        plot_kargs["xyLabels"] = ["$%s$" %(var),""]

    return plot_fcts(axis, x_vect, ys, **plot_kargs)


# It needs a grid with both x-offset and y-offset set to 1, resulting from ellipses_varsGrid()
def insert_gaussians_into_varsGrid(grid,matrices,vars=var_names,**plot_kargs):
    Nvars = len(vars)
    for i, var in enumerate(vars):
        if i==Nvars-1:
            plot_kargs["xyLabels"] = ["$%s$" %(var),""]
        else:
            plot_kargs["xyLabels"] = ["",""]
        ax_ii = pl.subplot(grid[i, i])
        gaussian_matplot(ax_ii,matrices,var,**plot_kargs)



def plot_fcts(axis, x, ys, **plot_kargs):
    plot_kargs = check_matplot_arguments("linePlot",**plot_kargs)

    matplotlib.rcParams['text.usetex'] = True
    font_style = {'weight' : 'normal', 'size': plot_kargs['ticks_size'],'family':'serif','serif':['Palatino']}
    matplotlib.rc('font',**font_style)

    for y, label, color, lineStyle, opacity in zip(ys,plot_kargs['labels'],plot_kargs['colors'],plot_kargs['lineStyles'],plot_kargs['opacities']):
        axis.plot(x,y,lineStyle+color,label=r'%s'%(label),alpha=opacity)
    if plot_kargs['grid']==True:
        axis.grid(True)
    axis.legend(loc=plot_kargs['legend'])
    if 'x' in plot_kargs['log']:
        if 'symx' in plot_kargs['log']:
            axis.set_xscale('symlog')
        else:
            axis.set_xscale('log')
    if 'symy' in plot_kargs['log']:
        axis.set_yscale('symlog')
    elif 'symx' in plot_kargs['log']:
        print "boh"
    elif 'y' in plot_kargs['log']:
        axis.set_yscale('log')
    if plot_kargs["xrange"]!=0:
        axis.set_xlim(plot_kargs["xrange"])
    if plot_kargs["yrange"]!=0:
        axis.set_ylim(plot_kargs["yrange"])
    axis.set_xlabel(r'%s' %(plot_kargs['xyLabels'][0]),fontsize=plot_kargs["label_size"])
    axis.set_ylabel(r'%s' %(plot_kargs['xyLabels'][1]),fontsize=plot_kargs["label_size"])

    return axis





# def plot_fcts(x, ys, pdf_name, labels, log="xy", legend="best",xyLabels=['$k$ [$h$/Mpc]', '$P(k)$ [(Mpc/$h$)$^3$]'],colors=['r','b','g'], lineStyles=['-']*15):
#     fig=pl.figure()
#     ax=fig.add_subplot(111)
#     for y, label, color, lineStyle in zip(ys,labels,colors,lineStyles):
#         ax.plot(x,y,lineStyle+color,label=r'%s'%(label))
#     ax.grid(True)
#     ax.legend(loc='best')
#     if 'x' in log:
#         if 'symx' in log:
#             ax.set_xscale('symlog')
#         else:
#             ax.set_xscale('log')
#     if 'y' in log:
#         if 'symy' in log:
#             ax.set_yscale('symlog')
#         else:
#             ax.set_yscale('log')
#     ax.set_xlabel(r'%s' %(xyLabels[0]),fontsize=15)
#     ax.set_ylabel(r'%s' %(xyLabels[1]),fontsize=15)
#     fig.set_tight_layout(True)
#     fig.savefig(pdf_name+".pdf")



def FOM(FM_matrices, vars_num):
    results = []
    for FM_matrix in FM_matrices:
        sub_matrix = np.linalg.inv(marginalise(FM_matrix,vars_num))
        results.append(np.sqrt(np.linalg.det(sub_matrix)))
    return results

def FOM_total(FM_matrices):
    results = []
    for FM_matrix in FM_matrices:
        # inv_matrix = np.linalg.inv(FM_matrix)
        results.append(np.sqrt(np.linalg.det(inv_matrix)))
    return results

