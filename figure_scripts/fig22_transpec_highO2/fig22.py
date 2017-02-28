"""
Figure 22:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig22"
    fname1 = "profile_o2lb_10bar_h2o.pt_filtered_transit.tran"
    fname2 = "profile_o2lb_10bar_dry.pt_filtered_transit.tran"
    title1 = "10 bar O$_2$ (Ocean)"
    title2 = "10 bar O$_2$ (Dessicated)"
    lammin = 0.2
    lammax = 20.
    ylim = [0.0, 100.0]
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
    fname = [os.path.join(os.path.dirname(__file__),"model_outputs/", f) for f in [fname1, fname2]]
    title = [title1, title2]
    labels = [labels1, labels2]

    # Make reflectance plot
    spc.plot_trans(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None)

    return
#########################################

if __name__ == "__main__":

    make_fig()
