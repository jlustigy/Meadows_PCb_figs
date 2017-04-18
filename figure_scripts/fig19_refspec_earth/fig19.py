"""
Figure 19:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig19"
    fname1 = "profile_Earth_proxb_.pt_hitran2012_o4_noh2co_187Kstrat_toa.rad"
    fname2 = "profile_Earth_proxb_.pt_stratocum_hitran2012_o4_noh2co_187Kstrat_toa.rad"
    fname3 = "profile_Earth_proxb_.pt_cirrus_hitran2012_o4_noh2co_187Kstrat_toa.rad"
    fname_new = "model_outputs/profile_Earth_proxb_combined_toa.rad"
    title1 = "Earth-like ($50\%$ clouds)"
    lammin = 0.1
    lammax = 2.5
    ylim = [-0.02, 0.3]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]

    moleloc = {"O$_2$" : [(0.15, 0.01), (0.63, 0.124), (0.68, 0.068), (0.77, 0.015), (1.275, .07)],
               "O$_3$" : [(0.25, 0.07), (0.62, 0.18)],
               "CO$_2$" : [(1.6, 0.12), (2.15, 0.09)],
               "CO" : [(0.15, 0.05)],
               "H$_2$O" : [(0.88, 0.04), (0.98, 0.19), (1.163, 0.16), (1.4, 0.12), (1.87, 0.08)],
               "CH$_4$" : [(1.22, 0.19), (1.4, 0.07), (1.65, 0.09), (1.87, 0.06), (2.05,0.06), (2.27,0.05), (2.37,0.04)]
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
