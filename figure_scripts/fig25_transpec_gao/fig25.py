"""
Figure 25:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig25"
    #fname = "Gao2015_case3.pt_filtered_transit.tran"
    fname = "smart_gao_1bar_FINAL.trnst"
    title = "1 bar CO$_2$/CO/O$_2$ (Desiccated)"
    lammin = 0.2
    lammax = 20.
    ylim = [0.0, 100.0]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]
    m = {
        "O$_2$" : [(0.69, 23.), (0.76, 25.), (0.88, 18.), (1.05, 21)],
        "O$_3$" : [(0.26, 72.),(0.6, 24.),(9.7,35.)],
        "CO" : [(1.6, 30.0), (2.35, 36.), (4.8, 46.)],
        #"CH$_4$" : [(0.9, 31.), (1.16,40.),(1.4, 43.), (2.33,49.), (3.4, 60), (6.6,39.), (7.9,53.)],
        #"H$_2$O" : [(0.86, 7.), (5.8, 34.)],
        "CO$_2$" : [(1.21, 28.), (2.01, 40.), (2.8, 48.), (4.3, 58.), (7.5, 25.), (15, 60.)]
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
                 title=title, labels=None, moleloc=m)

    return
#########################################

if __name__ == "__main__":

    make_fig()
