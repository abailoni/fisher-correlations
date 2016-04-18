#**************************************************************
#**************************************************************
#
# Import spectra from CLASS, normalise and take num. derivatives
#
#**************************************************************
#**************************************************************


PATH_SPECTRA = "INPUT/CLASS_spectra/"

k_extr = {}

def import_CLASS_data():
    print "\nInterpolating and normalising spectra. . .",
    global k_extr, norms, zero_spectrum_tools_anorm, zero_spectrum_tools, derivatives_tools_anorm, derivatives_tools, tools_anorm, tools
    for (dirpath, dirnames, filenames) in walk(PATH_SPECTRA):
        names = [ fi for fi in filenames if fi.endswith(".dat") ]
        break
    for filename in names:
        print".",
        # Determine the kind of data in the file:
        m = match(r"der_([^\.]+)", filename)
        if m:
            var = m.group(1)
        else: #this is the reference spectrum
            var = "spectrum"
        #
        #################################################
        # Adjust data because CLASS output file is a mess...
        # Import file line by line:
        with open(PATH_SPECTRA+filename) as fil:
            lines = fil.readlines()
        lines = lines[4:-1]
        vect_k = np.zeros(len(lines))
        data_class_anorm = np.zeros(len(lines))
        # Detect data in the line:
        for j in range(len(lines)):
            m = match(r" +([^ ]+) +([^ ]+)", lines[j])
            if m: vect_k[j], data_class_anorm[j] = float(m.group(1)), float(m.group(2))
        #
        #################################################
        # Interpolate: (GSL)
        if var=="spectrum":
            tools_anorm = &zero_spectrum_tools_anorm[0]
            tools = &zero_spectrum_tools[0]
        else: #var goes from 1 to 4
            tools_anorm = &derivatives_tools_anorm[n_var_import[var]]
            tools = &derivatives_tools[n_var_import[var]]
        alloc_interp_GSL(vect_k,data_class_anorm, tools_anorm)
        k_extr[var] = np.array([vect_k[1], vect_k[-2]])
        # Normalise:
        norms[var] = 2*np.pi**2 / quad( lambda k, var=var: k**3 * eval_interp_GSL(k, tools_anorm) * Fourier_W_k(k*R_sigma)**2 , k_extr[var][0], k_extr[var][1], full_output=1)[0]
        data_class_norm = data_class_anorm * sigma8_ref**2 * norms[var]
        alloc_interp_GSL(vect_k,data_class_norm, tools)
    # Adjust max and min for interpolation, because they differ in every spectrum:
    global k_min, k_max
    k_min, k_max = max([k_extr[var][0] for var in k_extr]), min([k_extr[var][1] for var in k_extr])
    print "Done!"

