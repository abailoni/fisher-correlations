import numpy as np
import matplotlib.pyplot as pl

def plot_fcts(x, ys, pdf_name, labels, log="xy", legend="best",xyLabels=['$k$ [$h$/Mpc]', '$P(k)$ [(Mpc/$h$)$^3$]'],colors=['r','b','g'], lineStyles=['-']*15):
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
santi_07 = np.loadtxt(open(INPUT_SANTI+"CAMB_spectrum_0.7.csv","rb"),delimiter=",",skiprows=0)
santi_0 = np.loadtxt(open(INPUT_SANTI+"CAMB_spectrum_0.csv","rb"),delimiter=",",skiprows=0)

from scipy.interpolate import interp1d
santi_0_interp = interp1d(santi_0[:,0],santi_0[:,1],kind="slinear")


# Half EUCLID-2012 bins:
N_bins = 7
N_k = 3000
k_vect = np.logspace(np.log10(1e-4),np.log10(1.4),N_k)
bin1, bin2 = 6, 6

cFM.init()
spectrum, santi = np.empty(N_k), np.empty(N_k)
for i in range(N_k):
    spectrum[i] = cFM.zero_spectrum_py(k_vect[i])
    santi[i] = santi_0_interp(k_vect[i])


# plot_fcts_PRO([santi_0[:,0], k_vect],[santi_0[:,1], spectrum], "difference_santi_0", ["Santi z=0", "Python wrapper"],'xy')
plot_fcts(k_vect,[(santi-spectrum)], "difference_santi_0", ["Difference: Santi-python"],'x')


