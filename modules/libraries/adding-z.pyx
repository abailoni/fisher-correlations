cdef:
    double[::1] Growth_der_redshift # N_bins



def compute_data_adding_z():
    global Growth_der_redshift
    Growth_der_redshift = Fun_der_redshift(Growth,z_avg)


def obsSpectr_redshift_der(bin1,bin2,k,mu,der_bin):
    return observed_spectrum(bin1,bin2,k,mu) * Growth_bins[der_bin] * Growth_der_redshift[der_bin]






