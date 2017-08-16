"""
Figure 01

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig01"
    fname = "fig35_modern_earth_standard_clearsky_2000_100000cm_60sza_sur.rad"
    title  = "Transmittance of the Modern Earth Atmosphere"

    lammin = 0.0
    lammax = 10.0
    ylim = [0.0, 1.0]

    # import basic dependencies
    import numpy as np
    import os
    import sys

    # import local package
    sys.path.append("../../")
    from utils import fig_params
    import utils.spectra as spc

    # More general params
    fname = os.path.join(os.path.dirname(__file__),"model_outputs/", fname)

    # Plot and save
    spc.plot_transmittance(fname, lammin=lammin, lammax=lammax, ylim=ylim,
                           title=title, savetag=savetag)

    """

    # Number of elements in a row
    Nrow = 14
    # Convert each line to vector, compose array of vectors
    arrays = np.array([np.array(map(float, line.split())) for line in open(fname)])
    # Flatten and reshape into rectangle grid
    arr = np.hstack(arrays).reshape((Nrow, -1), order='F')
    # Parse columns
    lam   = arr[0,:]
    wno   = arr[1,:]
    solar = arr[2,:]
    direct = arr[3,:]
    sky = arr[4,:]

    # Compute transmittance
    transmit = direct/solar


    import matplotlib as mpl
    import matplotlib.pyplot as plt
    fig_params.set_spectrum_figsize()

    fig, ax = plt.subplots()
    ax.plot(lam, transmit, lw=2.0, color="black")
    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel(r"Transmittance")
    ax.set_title(title)
    ax.set_xlim(lammin, lammax)
    ax.set_ylim(ylim)

    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')

    """

    return
#########################################

if __name__ == "__main__":

    make_fig()
