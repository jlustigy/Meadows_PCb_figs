"""
Figure 23:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig23"
    #fname1 = "old/profile_O2_CO2_10bar_prox_transit_hitran2012_50_100000cm.tran"
    #fname2 = "old/profile_O2_CO2_90bar_prox_transit_hitran2012_50_100000cm.tran"
    fname1 = "10bar_O2_CO2_final.pt_filtered_hitran2012_50_100000cm.trnst"
    fname2 = "90bar_O2_CO2_profile.pt_filtered_hitran2012_50_100000cm.trnst"
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
    m1 = {"O$_2$" : [(0.69, 38.), (0.76, 41.)],
        "O$_3$" : [(0.26, 71.),(0.6, 43.),(4.8,53.),(9.7,50.)],
            "O$_4$" : [(1.06, 34.), (1.27,38.)],
                "CO$_2$" : [(1.6, 41.), (2.0, 51.), (2.8, 59.), (4.3, 69.), (7.5, 39.), (15.5, 70.)]
                }
    m2 = {"O$_2$" : [(0.69, 57.5), (0.76, 61.)],
        "O$_3$" : [(0.26, 93.), (0.6, 60.),(4.8,73.),(9.7,70.)],
            "O$_4$" : [(1.06, 56.), (1.27,59.)],
            "CO$_2$" : [(1.6, 63.), (2.0, 71.), (2.8, 80.), (4.3, 89.), (7.5, 60.), (15.5, 89.)]
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
    """
    spc.plot_trans(fname, savetag=savetag+"2", lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, moleloc=moleloc)
    """

    title = [r"O$_2$-CO$_2$", ""]

    spc.plot_trans(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, forced_single=True,
                 plot_kwargs=plot_kwargs, legend=True, moleloc=[None, m2])

    return
#########################################

if __name__ == "__main__":

    make_fig()
