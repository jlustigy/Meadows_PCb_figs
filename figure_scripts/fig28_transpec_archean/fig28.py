"""
Figure 28:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig28"
    fname1 = "HAZE_1.50e-02ch4_clear_new.trnst"
    fname2 = "HAZE_1.00e-02ch4_clear_new.trnst"
    title  = "Archean Earth-like"
    title1 = "Haze, $5\%$ CO$_2$, $1.5\%$ CH$_4$"
    title2 = "No Haze, $5\%$ CO$_2$, $1\%$ CH$_4$"
    plot_kwargs = [
        {"color" : "orange", "label" : title1, "alpha" : 0.7},
        {"color" : "purple", "label" : title2, "alpha" : 0.7}
    ]
    lammin = 0.21
    lammax = 20.
    ylim = [0.0, 100.0]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]

    moleloc = {
                "CH$_4$" : [(0.73, 29.), (0.79, 7.0), (0.9, 48.), (1.0, 37.), (1.16,55.), (1.4, 60.), (2.33,61.), (3.4, 73), (6.6,51.), (7.9,63.)],
                "CO$_2$" : [(2.0, 47.), (2.8, 56.), (4.3, 70.), (9.5, 32.), (10.4, 27.), (15., 73.)],
                "CO" : [(2.6, 20.), (4.8, 47.)],
                "C$_2$H$_6$" : [(7.0, 20.), (12., 38.)],
                "Haze" : [(0.4, 63.), (5.85, 45.5)]
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
    fname = [os.path.join(os.path.dirname(__file__),"model_outputs/", f) for f in [fname1, fname2]]
    #title = [title1, title2]
    title = [title, ""]
    labels = [labels1, labels2]

    # Make reflectance plot
    spc.plot_trans(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, forced_single=True,
                 plot_kwargs=plot_kwargs, legend=True, moleloc=[moleloc,None])

    return
#########################################

if __name__ == "__main__":

    make_fig()
