"""
Figure 26:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig26"
    #fname = "profile_earth_prox.pt_filtered_transit_hitran2012_50_100000cm.tran"
    #fname = "profile_Earth_proxb_clear.trnst"
    fname = "profile_Earth_proxb_.pt_hitran2012_o4_noh2co_187Kstrat.trnst"
    title = "Earth-like"
    lammin = 0.2
    lammax = 20.
    ylim = [0.0, 100.0]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]

    moleloc = {
                "O$_2$" : [(0.215, 57.), (0.69, 24.), (0.76, 27.)],
                "O$_3$" : [(0.28, 60.),(0.6, 25.),(4.8,43.),(9.7,39.)],
                "CH$_4$" : [(0.9, 31.), (1.0, 26.), (1.16,40.),(1.4, 43.), (2.33,49.), (3.4, 60), (6.6,39.), (7.9,53.)],
                "H$_2$O" : [(0.86, 7.), (5.8, 34.)],
                "CO$_2$" : [(1.92, 40.), (2.8, 51.), (4.3, 61.),(15.5, 65.)]
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
    fname = os.path.join(os.path.dirname(__file__),"model_outputs/", fname)

    # Make reflectance plot
    spc.plot_trans(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, moleloc=moleloc)

    return
#########################################

if __name__ == "__main__":

    make_fig()
