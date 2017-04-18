"""
Figure 18:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig18"
    #fname1 = "gao_1bar_final_toa_old.rad"
    fname1 = "smart_gao_1bar_FINAL_toa.rad"
    title1 = "1 bar CO$_2$/CO/O$_2$ (Dessicated)"
    lammin = 0.1
    lammax = 2.5
    ylim = [-0.01, 0.3]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]


    moleloc = {"O$_2$" : [(0.15, 0.01), (0.63, 0.10), (0.69, 0.05), (0.76, -0.001), (1.27, .057)],
               "O$_3$" : [(0.25, 0.01), (0.55, 0.125)],
               "CO$_2$" : [(0.88, 0.065), (1.05, 0.022), (1.21, .175), (1.62, 0.13), (1.75, 0.17), (2.01, 0.13), (2.45, 0.15)],
               "CO" : [(0.15, 0.03), (1.2, 0.19), (1.58, 0.17), (2.35, 0.13)]
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
    fname = os.path.join(os.path.dirname(__file__),"model_outputs/", fname1)

    # Make reflectance plot
    spc.plot_rad(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title1, labels=None, moleloc=moleloc)

    return
#########################################

if __name__ == "__main__":

    make_fig()
