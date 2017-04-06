"""
Figure 22:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig22_new"
    #fname1 = "profile_o2lb_10bar_h2o.pt_filtered_transit.tran"
    #fname2 = "profile_o2lb_10bar_dry.pt_filtered_transit.tran"
    fname1 = "10bar_O2_wet.trnst"
    fname2 = "10bar_O2_dry.trnst"
    title1 = "10 bar O$_2$ (Ocean)"
    title2 = "10 bar O$_2$ (Dessicated)"
    lammin = 0.2
    lammax = 20.
    ylim = [0.0, 100.0]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]
    plot_kwargs = [
        {"color" : "purple", "label" : "Ocean", "alpha" : 0.7},
        {"color" : "orange", "label" : "Dessicated", "alpha" : 0.7}
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

    title = [r"10 bar High O$_2$", ""]

    spc.plot_trans(fname, savetag=savetag+"1", lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, forced_single=True,
                 plot_kwargs=plot_kwargs, legend=True)

    return
#########################################

if __name__ == "__main__":

    make_fig()
