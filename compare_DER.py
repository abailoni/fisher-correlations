import numpy as np
import matplotlib.pyplot as pl

def plot_fcts(x, ys, pdf_path, labels, log="xy", title="", legend="best",xyLabels=['$k$ [$h$/Mpc]', '$P(k)$ [(Mpc/$h$)$^3$]'],colors=['r','b','g','y','p'], lineStyles=['-']*15):
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
    ax.set_title(title)
    ax.set_ylabel(xyLabels[1])
    fig.set_tight_layout(True)
    fig.savefig(pdf_path+".pdf")

def plot_fcts_PRO(xs, ys, pdf_path, labels, log="xy", title="", legend="best",xyLabels=['$k$ [$h$/Mpc]', '$P(k)$ [(Mpc/$h$)$^3$]'],colors=['r','b','g'], lineStyles=['-']*15):
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
    ax.set_title(title)
    ax.set_xlabel(xyLabels[0])
    ax.set_ylabel(xyLabels[1])
    fig.set_tight_layout(True)
    fig.savefig(pdf_path+".pdf")

INPUT_SANTI = "santiDER/"
OWN_DIR = "ownDER/"

var_names = ['$\\Omega_\\mathrm{b}$', '$\\Omega_\\mathrm{c}$', '$w_0$']
quantities = ["Hubble", "Growth", "D_A", "Beta"]
santiDER = [ np.loadtxt(open(INPUT_SANTI+name+".txt","rb"),delimiter=";") for name in quantities]
import modules.cFM as cFM
cFM.init()
ownDER = cFM.print_DER()
# ownDER = santiDER

# # Beta log:
# for bin in range(14):
#     ownDER[3][:,bin] = ownDER[3][:,bin]*cFM.beta(bin)

bins = np.arange(14)
perc = [ (santiDER[i]-ownDER[i])/santiDER[i]*100. for i in range(len(quantities))]

# print perc[0].shape


# plot_fcts_PRO([santi_0[:,0], k_vect],[santi_0[:,1], spectrum], "difference_santi_0", ["Santi z=0", "Python wrapper"],'xy')

for i, quantity in enumerate(quantities):
    plot_fcts(bins,perc[i], "plot_comparison_"+quantity, var_names,'', "Quantity: "+quantity, 'best',  ['Bins', 'Percentage difference [wrt Santi]'])


