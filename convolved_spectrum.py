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
        ax.set_yscale('log')
    ax.set_xlabel(r'%s' %(xyLabels[0]),fontsize=15)
    ax.set_ylabel(r'%s' %(xyLabels[1]),fontsize=15)
    fig.set_tight_layout(True)
    fig.savefig(OUTPUT_DIR+pdf_name+".pdf")

OUTPUT_DIR = "plots/convolved_spectra/"

import modules.cFM as cFM

# Half EUCLID-2012 bins:
N_bins = 7
N_k = 3000
k_vect = np.logspace(np.log10(1e-4),np.log10(1.),N_k)
bin1, bin2 = 0, 0

cFM.init()
spectrum, convSpec_14 = np.empty(N_k), np.empty(N_k)
for i in range(N_k):
    convSpec_14[i] = cFM.growth_spectrum_py(bin1,bin2,k_vect[i],"windowFun")
    spectrum[i] = cFM.growth_spectrum_py(bin1,bin2,k_vect[i],"uncorr")

z_in = [0.65, 0.85, 1.05, 1.25, 1.45, 1.65, 1.85, 2.05]
z_avg = [ (z_in[i]+z_in[i+1])/2. for i in range(N_bins)]
dens = cFM.density_py(np.array(z_avg))
bias = cFM.bias_py(z_avg)
cFM.set_survey(bins_list=z_in, dens_list=dens,bias_list=bias)
cFM.init()
convSpec_7 = np.empty(N_k)
for i in range(N_k):
    convSpec_7[i] = cFM.growth_spectrum_py(bin1,bin2,k_vect[i],"windowFun")


plot_fcts(k_vect[:-20],[spectrum[:-20],convSpec_7[:-20],convSpec_14[:-20]], "compare_bins_DEF_%d_%d" %(bin1,bin2), ["$P_{\\mathrm{lin}}(k;z_{%d})$" %(bin1+1), "$P_{\\mathrm{conv}}(k;z_{%d})$ with bin width $\\Delta z=2$" %(bin1+1), "$P_{\\mathrm{conv}}(k;z_{%d})$ with bin width $\\Delta z=1$" %(bin1+1)],"xy", "best", ['$k$ [$h$/Mpc]', '$P(k)$ [(Mpc/$h$)$^3$]'], ['r','g','b'])

# plot_fcts(k_vect,[convSpec_14 - convSpec_7], "difference_%d_%d" %(bin1,bin2), ["$P_{%d,%d}^{14}-P_{%d,%d}^{7}$" %(bin1,bin2,bin1,bin2)],'x')


