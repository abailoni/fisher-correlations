import numpy as np
import matplotlib
import matplotlib.pyplot as pl

def print_matrix_with_values(matrix, output_pdf_path='matrix_plot.pdf', plot_title='',font_size='8', font_size_title='10', cmap='jet'):
    font_ticks = {'size': font_size}
    matplotlib.rc('font', **font_ticks)

    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.set_title(plot_title, fontsize=font_size_title)

    # Plot color matrix:
    cax = ax.matshow(matrix,cmap=cmap)
    cbar = fig.colorbar(cax)
    ax.set_aspect('equal')

    # Text portion:
    ind_x, ind_y = np.arange(matrix.shape[0]), np.arange(matrix.shape[1])
    x, y = np.meshgrid(ind_x, ind_y)
    for x_val, y_val in zip(x.flatten(), y.flatten()):
        ax.text(x_val, y_val, "%.1f" %(matrix[x_val,y_val]), va='center', ha='center')

    fig.savefig(output_pdf_path)



PATH_IN = "OUTPUT/FFT/"
PATH_PDF = "plots/"
file_name = "FMcorr_EUCLID2012_AP1-14bins"

matrix = np.loadtxt(open(PATH_IN+file_name+".csv"),delimiter=" ",skiprows=0)

# --- Convert to symlog ---
# idxs_pos, idxs_neg = (matrix>0).nonzero(), (matrix<0).nonzero()
# matrix[idxs_pos], matrix[idxs_neg] = np.log10(matrix[idxs_pos]), -np.log10(-matrix[idxs_neg])

print_matrix_with_values(matrix, PATH_PDF+file_name+".pdf","EUCLID-2012, with AP, 14 bins", '3', '11')
