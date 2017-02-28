"""
Figure 19:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig19"
    fname1 = "profile_earth_prox.pt_filtered_hitran2012_50_100000cm_toa.rad"
    fname2 = "profile_earth_prox.pt_stratocum_hitran2012_50_100000cm_toa.rad"
    fname3 = "profile_earth_prox.pt_cirrus_hitran2012_50_100000cm_toa.rad"
    fname_new = "model_outputs/profile_earth_prox.pt_combined_hitran2012_50_100000cm_toa.rad"
    title1 = "Earth-like ($50\%$ clouds)"
    lammin = 0.1
    lammax = 2.5
    ylim = [-0.02, 0.3]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]

    # import basic dependencies
    import numpy as np
    import os
    import sys

    # import local package
    sys.path.append("../../")
    from utils import fig_params
    import utils.spectra as spc

    # More general params
    fname = [os.path.join(os.path.dirname(__file__),"model_outputs/", f) for f in [fname1, fname2, fname3]]

    # Weight each flux
    dic = {
        fname[0] : 0.5,
        fname[1] : 0.25,
        fname[2] : 0.25
    }

    # Weight spectra and create output file
    spc.weight_spectra(dic, output=fname_new, ret=False)

    # Make reflectance plot
    spc.plot_rad(fname_new, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                     title=title1, labels=None)

    return
#########################################

if __name__ == "__main__":

    make_fig()
