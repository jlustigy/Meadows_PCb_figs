"""
Figure 16:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig16"
    fname1 = "profile_O2_CO2_10bar_prox_hitran2012_50_100000cm_toa.rad"
    fname2 = "profile_O2_CO2_90bar_prox_hitran2012_50_100000cm_toa.rad"
    title1 = "10 bar O$_2$-CO$_2$"
    title2 = "90 bar O$_2$-CO$_2$"
    lammin = 0.2
    lammax = 2.5
    ylim = [-0.02, 0.51]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]

    m1 = {"O$_2$" : [(0.23, 0.1), (0.67, 0.21), (0.76, 0.21)],
          "CO" : [(0.25, 0.05), (2.35, 0.06)],
          "O$_3$" : [(0.28, 0.15)],
          "O$_4$" : [(0.39, 0.39), (0.46, 0.33), (0.49, 0.29), (0.55, 0.2), (0.58, 0.0), (0.64, 0.0), (1.06, 0.1), (1.27, 0.1)],
          "CO$_2$" : [(0.88, 0.23), (1.6, 0.08), (1.8, .23), (2.01, 0.025), (2.2, 0.13)]
          }
    m2 = {"O$_2$" : [(0.23, 0.12), (0.67, 0.23), (0.76, 0.23)],
          "CO" : [(0.25, 0.05), (2.35, 0.025)],
          "O$_3$" : [(0.28, 0.15)],
          "O$_4$" : [(0.39, 0.48), (0.47, 0.42), (0.49, 0.35), (0.55, 0.2), (0.58, 0.0), (0.64, 0.0), (1.06, 0.025), (1.27, 0.025)],
          "CO$_2$" : [(0.88, 0.32), (1.6, 0.025), (1.8, .025), (2.01, 0.025), (2.2, 0.025)]
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
