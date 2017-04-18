"""
Figure 16:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig16"
    fname1 = "10bar_O2_CO2_final.pt_filtered_hitran2012_50_100000cm_toa.rad"
    fname2 = "90bar_O2_CO2_profile.pt_filtered_hitran2012_50_100000cm_toa.rad"
    title1 = "10 bar O$_2$-CO$_2$"
    title2 = "90 bar O$_2$-CO$_2$"
    lammin = 0.2
    lammax = 2.5
    ylim = [-0.02, 0.51]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]

    m1 = {"O$_2$" : [(0.23, 0.1), (0.69, 0.25), (0.76, 0.25)],
          "CO" : [(0.25, 0.05), (2.35, 0.11)],
          "O$_3$" : [(0.28, 0.15)],
          "O$_4$" : [(0.39, 0.44), (0.46, 0.37), (0.49, 0.34), (0.55, 0.25), (0.58, 0.03), (0.64, 0.01), (1.06, 0.08), (1.27, 0.1)],
          "CO$_2$" : [(0.88, 0.25), (1.6, 0.08), (1.8, .23), (2.01, 0.025), (2.2, 0.13)]
          }
    m2 = {"O$_2$" : [(0.23, 0.12), (0.69, 0.36), (0.762, 0.26)],
          "CO" : [(0.25, 0.05), (2.35, 0.11)],
          "O$_3$" : [(0.28, 0.15)],
          "O$_4$" : [(0.39, 0.495), (0.47, 0.47), (0.49, 0.165), (0.55, 0.36), (0.58, 0.03), (0.64, 0.01), (1.06, 0.025), (1.27, 0.025)],
          "CO$_2$" : [(0.88, 0.355), (1.6, 0.025), (1.8, .025), (2.01, 0.04), (2.2, 0.07)]
          }
    moleloc = [m1, m2]

    # import basic dependencies
    import numpy as np
    import os
    import sys

    # import local package
    sys.path.append("../../")
    from utils import fig_params
    import utils.spectra as spc

    # More general params
    fname = [os.path.join(os.path.dirname(__file__),"model_outputs/", f) for f in [fname1, fname2]]
    title = [title1, title2]
    labels = [labels1, labels2]

    # Make reflectance plot
    spc.plot_rad(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, moleloc=moleloc)

    return
#########################################

if __name__ == "__main__":

    make_fig()
