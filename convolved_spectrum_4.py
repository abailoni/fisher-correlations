import numpy as np
import matplotlib.pyplot as pl
import matplotlib

matplotlib.rcParams['text.usetex'] = True
# matplotlib.rc('font',**{'family':'serif'})
matplotlib.rc('font',**{'family':'serif','serif':['Palatino'],'size': 15})

def plot_fcts(x, ys, pdf_name, labels, log="xy", legend="best",xyLabels=['$k$ [$h$/Mpc]', '$P(k)$ [(Mpc/$h$)$^3$]'],colors=['r','b','g'], lineStyles=['-']*15):
    fig=pl.figure()
    ax=fig.add_subplot(111)
    for y, label, color, lineStyle in zip(ys,labels,colors,lineStyles):
        ax.plot(x,y,lineStyle+color,label=r'%s'%(label))
    ax.grid(True)
    ax.legend(loc='best')
    if 'x' in log:
        ax.set_xscale('log')
    if 'y' in log:
        ax.set_yscale('symlog')
    ax.set_xlabel(r'%s' %(xyLabels[0]),fontsize=15)
    ax.set_ylabel(r'%s' %(xyLabels[1]),fontsize=15)
    fig.set_tight_layout(True)
    fig.savefig(OUTPUT_DIR+pdf_name+".pdf")

OUTPUT_DIR = "plots/convolved_spectra/"

import modules.cFM as cFM

N_k = 3000
k_vect = np.logspace(np.log10(1e-4),np.log10(1.),N_k)
bin1, bin2 = 0, 13

cFM.init()
spectrum, convSpec_1, convSpec_2, convSpec_3 = np.empty(N_k), np.empty(N_k), np.empty(N_k), np.empty(N_k)
for i in range(N_k):
    convSpec_1[i] = cFM.growth_spectrum_py(0,1,k_vect[i],"windowFun")
    convSpec_2[i] = cFM.growth_spectrum_py(0,0,k_vect[i],"windowFun")
    # convSpec_3[i] = cFM.growth_spectrum_py(0,13,k_vect[i],"windowFun")
    spectrum[i] = cFM.growth_spectrum_py(0,0,k_vect[i],"uncorr")


plot_fcts(k_vect[:-80],[spectrum[:-80], convSpec_2[:-80],convSpec_1[:-80]], "compare_corr_01", ["$P_{\\mathrm{lin}}(k;z_1)$","$P_{\\mathrm{conv}}(k;z_1)$","$P_{\\mathrm{conv}}(k;z_1,z_2)$"])

# plot_fcts(k_vect,[convSpec_14 - convSpec_7], "difference_%d_%d" %(bin1,bin2), ["$P_{%d,%d}^{14}-P_{%d,%d}^{7}$" %(bin1,bin2,bin1,bin2)],'x')


