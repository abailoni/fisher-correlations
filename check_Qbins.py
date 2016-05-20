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


OUTPUT_DIR = "plots/convolved_spectra/"

import modules.cFM as cFM

N_k = 3000
k_vect = np.logspace(np.log10(1e-4),np.log10(1.),N_k)
bin1, bin2 = 0, 0

cFM.compute_CAMB_spectra()
cFM.import_zero_spectrum_der_k()
cFM.compute_survey_DATA()
cFM.store_int1_test()
cFM.set_typeFM("windowFun")

spectrum, convSpect, K = np.empty(N_k), np.empty(N_k), np.empty(N_k)
for i in range(N_k):
    if k_vect[i]>0.001 and k_vect[i]<0.1:
        spectrum[i] = 1.
    else:
        spectrum[i] = 0.
    K[i] = cFM.K_py(k_vect[i],bin1,bin2)
    convSpect[i]= cFM.spectrum(k,bin1,bin2)

k1, k2 = np.linspace(0.00000001,1.4,7000), np.linspace(0.00000001,10,70000)
convSpect1 = cFM.test_integral_1(bin1,bin2,"spectrum",k1)
convSpect2 = cFM.test_integral_1(bin1,bin2,"spectrum",k2)

np.savetxt("Q_constant_01.csv",np.column_stack((k_vect,K)))


k_deciso = 0.05
idx1, idx2 = (np.abs(k1-k_deciso)<5e-4).nonzero(), (np.abs(k2-k_deciso)<5e-4).nonzero()
print convSpect1[idx1[0][0]],convSpect2[idx2[0][0]]
print (convSpect1[idx1[0][0]]-convSpect2[idx2[0][0]])/convSpect1[idx1[0][0]]*100

plot_fcts_PRO([k1,k2],[convSpect1,convSpect2],"test_again",["1","2"],"x")
quit()

np.savetxt("K_01.csv",np.column_stack((k_vect,K)))

plot_fcts(k_vect[:-80],[spectrum[:-80],convSpect[:-80], K[:-80]], "test_Q_bins%d_%d" %(bin1,bin2), ["$P_{\\mathrm{lin}}(k;z_{%d})$" %(bin1+1), "$P_{\\mathrm{conv}}(k;z_{%d})$" %(bin1+1), "Qij"],"x")

# plot_fcts(k_vect,[convSpec_14 - convSpec_7], "difference_%d_%d" %(bin1,bin2), ["$P_{%d,%d}^{14}-P_{%d,%d}^{7}$" %(bin1,bin2,bin1,bin2)],'x')


