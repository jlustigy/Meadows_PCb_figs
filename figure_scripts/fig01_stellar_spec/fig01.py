"""
Figure 1

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    import numpy as np
    import os
    import sys
    sys.path.append("../../")
    from utils import fig_params
    import matplotlib.pyplot as plt
    from scipy.signal import savgol_filter

    # Params specific to this plot
    savetag = "fig01"
    plotdir="../../figures/"
    fname1 = "proxima_cen_1AU.txt"
    fname2 = "Kurucz1cm-1_susim_atlas2.dat"

    # Smoothing parameters
    window1 = 11
    window2 = 11
    poly1 = 5
    poly2 = 5

    # Read in atm files
    path = os.path.join(os.path.dirname(__file__),"data_files/", fname1)
    data = np.genfromtxt(path, skip_header=37)
    wl_prox, flux_prox = data[:,0], data[:,1]*(1.0/0.0485)**2
    #flux_prox = savgol_filter(flux_prox, window1, poly1, mode='nearest')

    path = os.path.join(os.path.dirname(__file__),"data_files/", fname2)
    data = np.genfromtxt(path, skip_header=12)
    wl_sun, flux_sun = 1e4/data[:,0], data[:,1]*1e-4*data[:,0]**2
    #flux_sun = savgol_filter(flux_sun, window2, poly2, mode='nearest')

    fig, ax = plt.subplots(2, figsize=(10,10))

    ax[0].set_ylabel("Flux [W/m$^2$/$\mu$m]")
    ax[0].set_xlabel("Wavelength [$\mu$m]")
    ax[0].plot(wl_prox, flux_prox, color="orange", label="Proxima Centauri \nat 0.0485 AU", zorder=10)
    ax[0].plot(wl_sun, flux_sun, color="black", label="Sun at 1 AU")
    ax[0].set_xlim([0.0, 2.6])
    ax[0].set_ylim([0,2600])

    ax[1].set_ylabel("Flux [W/m$^2$/$\mu$m]")
    ax[1].set_xlabel("Wavelength [$\mu$m]")
    ax[1].semilogy()
    ax[1].plot(wl_prox, flux_prox, color="orange", label="Proxima Centauri at 0.0485 AU")
    ax[1].plot(wl_sun, flux_sun, color="black", label="Sun at 1 AU")
    ax[1].set_xlim([0.1, 0.4])
    #ax[1].set_ylim([0,2600])

    # Label
    leg=ax[0].legend(loc=0, fontsize=16)
    leg.get_frame().set_alpha(0.0)

    # Save figure
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')

    return
#########################################

if __name__ == "__main__":

    make_fig()