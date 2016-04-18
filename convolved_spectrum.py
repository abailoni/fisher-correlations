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
    fig.savefig(OUTPUT_DIR+pdf_name+".pdf")

OUTPUT_DIR = "plots/convolved_spectra/"

import modules.cFM as cFM

# Half EUCLID-2012 bins:
N_bins = 7
N_k = 3000
k_vect = np.linspace(1e-4,1.4,N_k)
bin1, bin2 = 6, 6

cFM.init()
spectrum, convSpec_14 = np.empty(N_k), np.empty(N_k)
for i in range(N_k):
    convSpec_14[i] = cFM.windowed_Spectrum_py(k_vect[i],bin1,bin2)
    spectrum[i] = cFM.zero_spectrum_py(k_vect[i])

z_in = [0.65, 0.85, 1.05, 1.25, 1.45, 1.65, 1.85, 2.05]
z_avg = [ (z_in[i]+z_in[i+1])/2. for i in range(N_bins)]
dens = cFM.density_py(np.array(z_avg))
bias = cFM.bias_py(z_avg)
cFM.set_survey(bins_list=z_in, dens_list=dens,bias_list=bias)
cFM.init()
convSpec_7 = np.empty(N_k)
for i in range(N_k):
    convSpec_7[i] = cFM.windowed_Spectrum_py(k_vect[i],bin1,bin2)

# plot_fcts(k_vect,[spectrum,convSpec_7,convSpec_14], "compare_bins_%d_%d" %(bin1,bin2), ["CAMB spectrum", "7 bins $P_{%d,%d}$" %(bin1,bin2), "14 bins $P_{%d,%d}$" %(bin1,bin2)])

plot_fcts(k_vect,[convSpec_14 - convSpec_7], "difference_%d_%d" %(bin1,bin2), ["$P_{%d,%d}^{14}-P_{%d,%d}^{7}$" %(bin1,bin2,bin1,bin2)],'x')


