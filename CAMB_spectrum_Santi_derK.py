import numpy as np
import matplotlib.pyplot as pl

def plot_fcts(x, ys, pdf_name, labels, log="xy", legend="best",xyLabels=['$k$ [$h$/Mpc]', '$P(k)$ [(Mpc/$h$)$^3$]'], ylim=0, colors=['r','b','g'], lineStyles=['-']*15):
    fig=pl.figure()
    ax=fig.add_subplot(111)
    for y, label, color, lineStyle in zip(ys,labels,colors,lineStyles):
        ax.plot(x,y,lineStyle+color,label=label)
    ax.grid(True)
    ax.legend(loc='best')
    if 'x' in log:
        ax.set_xscale('log')
    if 'y' in log:
        ax.set_yscale('log')
    if ylim!=0:
        ax.set_ylim(ylim)
    ax.set_xlabel(xyLabels[0])
    ax.set_ylabel(xyLabels[1])
    fig.set_tight_layout(True)
    fig.savefig(OUTPUT_DIR+pdf_name+".pdf")

def plot_fcts_PRO(xs, ys, pdf_name, labels, log="xy", legend="best",xyLabels=['$k$ [$h$/Mpc]', '$P(k)$ [(Mpc/$h$)$^3$]'],colors=['r','b','g'], lineStyles=['-']*15):
    fig=pl.figure()
    ax=fig.add_subplot(111)
    for x, y, label, color, lineStyle in zip(xs,ys,labels,colors,lineStyles):
        ax.plot(x,y,lineStyle+color,label=label)
    ax.grid(True)
    ax.legend(loc='best')
    if 'x' in log:
        ax.set_xscale('log')
    if 'y' in log:
        ax.set_yscale('log')
    ax.set_xlabel(xyLabels[0])
    ax.set_ylabel(xyLabels[1])
    fig.set_tight_layout(True)
    fig.savefig(OUTPUT_DIR+pdf_name+".pdf")


INPUT_SANTI = "INPUT/Santi/"
OUTPUT_DIR = "plots/convolved_spectra/"

import modules.cFM as cFM


# Santi power spectra:
# santi_07 = np.loadtxt(open(INPUT_SANTI+"CAMB_spectrum_0.7.csv","rb"),delimiter=",",skiprows=0)
santi_0 = np.loadtxt(open(INPUT_SANTI+"CAMB_derK_0.csv","rb"),delimiter="\t",skiprows=0)

from scipy.interpolate import interp1d
from scipy.misc import derivative
santi_0_interp = interp1d(np.exp(santi_0[:,0]),santi_0[:,1],kind="slinear")


N_k = 3000
k_vect = np.logspace(np.log10(1e-4),np.log10(1.4),N_k)

cFM.init()
own_derK, santi = np.empty(N_k), santi_0_interp(k_vect)
for i in range(N_k):
    own_derK[i] = cFM.zero_spectrum_der_k_py(k_vect[i])/cFM.zero_spectrum_py(k_vect[i])
# own_interp = interp1d(k_vect,spectrum,kind="linear")

# k_vect_DER = k_vect[50:-50]
# own_derK = (lambda k: derivative(own_interp,k,1e-5))(k_vect_DER)
# santi_derK = (lambda k: derivative(santi_0_interp,k,1e-5))(k_vect_DER)
perc = (santi-own_derK)/santi*100

plot_fcts(k_vect,[own_derK,santi], "comparison_plots_derK", ["Python", "Cosmomathica"], 'xy', 'best', ['$k$', '$dP/dk$'])

# plot_fcts_PRO([santi_0[:,0], k_vect],[santi_0[:,1], spectrum], "difference_santi_0", ["Santi z=0", "Python wrapper"],'xy')
plot_fcts(k_vect,[(santi-own_derK)], "comparison_derK_abs", ["Absolute difference derivative k: Santi-python"],'x', 'best', ['$k$', 'Absolute difference $dP/dk$'])
plot_fcts(k_vect,[perc], "comparison_derK_der", ["Relative difference derivative k: Santi-python"],'x', 'best', ['$k$', 'Relative %% difference $dP/dk$'], [-30,30], ['g','b','g'], ['.'])

# plot_fcts(k_vect_DER,[(santi_derK-own_derK)/santi_derK*100], "comparison_santi-spectrum-perc_z0_derK", ["Difference derivative k: Santi-python (percentage)"],'x', 'best', ['$k$', '$dP/dk$'] )

# quit()

# --- Convert to symlog ---
# idxs_pos, idxs_neg = (perc>0).nonzero(), (perc<0).nonzero()
# perc[idxs_pos], perc[idxs_neg] = np.log10(perc[idxs_pos]), -np.log10(-perc[idxs_neg])

# plot_fcts(k_vect,[np.log(np.abs(perc))], "comparison_santi-spectrum-perc_z0_derK_symlog", ["Percentage difference dP/dk [wrt cosmomathica] - symLog "],'x', 'best', ['$k$', 'Log percentage difference Log$_{10}(dP/dk)$'])


