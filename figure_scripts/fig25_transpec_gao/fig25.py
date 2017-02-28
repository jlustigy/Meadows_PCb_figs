"""
Figure 25:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig25"
    #fname = "Gao2015_case3.pt_filtered_transit.tran"
    fname = "smart_gao_1bar_update_xsec.trnst"
    title = "1 bar CO$_2$/CO/O$_2$ (Dessicated)"
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
    fname = os.path.join(os.path.dirname(__file__),"model_outputs/", fname)

    # Make reflectance plot
    spc.plot_trans(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None)

    return
#########################################

if __name__ == "__main__":

    make_fig()
