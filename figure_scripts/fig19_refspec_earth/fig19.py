"""
Figure 19:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig19_new"
    fname1 = "profile_Earth_proxb_clear_toa.rad"
    fname2 = "profile_Earth_proxb_stratocum_toa.rad"
    fname3 = "profile_Earth_proxb_cirrus_toa.rad"
    fname_new = "model_outputs/profile_Earth_proxb_combined_toa.rad"
    title1 = "Earth-like ($50\%$ clouds)"
    lammin = 0.1
    lammax = 2.5
    ylim = [-0.02, 0.3]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]

    moleloc = {"O$_2$" : [(0.15, 0.01), (0.63, 0.04), (0.69, 0.019), (0.76, 0.0), (1.27, .013)],
               "O$_3$" : [(0.25, 0.01), (0.62, 0.065)],
               "CO$_2$" : [(0.88, 0.05), (1.05, 0.045), (1.21, .046), (1.62, 0.03), (1.75, 0.051), (2.01, 0.03), (2.45, 0.055)],
               "CO" : [(0.15, 0.02), (1.2, 0.054), (1.58, 0.047), (2.35, 0.03)]
          }

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
                     title=title1, labels=None, moleloc=moleloc)

    return
#########################################

if __name__ == "__main__":

    make_fig()
