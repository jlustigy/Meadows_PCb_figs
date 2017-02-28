"""
Figure 18:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig18"
    #fname1 = "gao_1bar_final_toa_old.rad"
    fname1 = "smart_gao_1bar_update_xsec_toa.rad"
    title1 = "1 bar CO$_2$/CO/O$_2$ (Dessicated)"
    lammin = 0.1
    lammax = 2.5
    ylim = [-0.005, 0.12]
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
    fname = os.path.join(os.path.dirname(__file__),"model_outputs/", fname1)

    # Make reflectance plot
    spc.plot_rad(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title1, labels=None)

    return
#########################################

if __name__ == "__main__":

    make_fig()
