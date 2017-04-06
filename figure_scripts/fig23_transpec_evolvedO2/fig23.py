"""
Figure 23:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig23_new"
    #fname1 = "old/profile_O2_CO2_10bar_prox_transit_hitran2012_50_100000cm.tran"
    #fname2 = "old/profile_O2_CO2_90bar_prox_transit_hitran2012_50_100000cm.tran"
    fname1 = "10bar_O2_CO2.trnst"
    fname2 = "90bar_O2_CO2.trnst"
    title1 = "10 bar O$_2$-CO$_2$"
    title2 = "90 bar O$_2$-CO$_2$"

    lammin = 0.2
    lammax = 20.
    ylim = [0.0, 100.0]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]
    plot_kwargs = [
        {"color" : "purple", "label" : "10 bar", "alpha" : 0.7},
        {"color" : "orange", "label" : "90 bar", "alpha" : 0.7}
    ]

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
    spc.plot_trans(fname, savetag=savetag+"2", lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None)

    title = [r"O$_2$-CO$_2$", ""]

    spc.plot_trans(fname, savetag=savetag+"1", lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, forced_single=True,
                 plot_kwargs=plot_kwargs, legend=True)

    return
#########################################

if __name__ == "__main__":

    make_fig()
