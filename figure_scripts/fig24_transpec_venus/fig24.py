"""
Figure 24:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig24"
    #fname1 = "fig24_tran_smart_spectra_pandora10bar_cloudy_500_100000cm-1.tran"
    #fname2 = "fig24_tran_smart_spectra_pandora90bar_clouds_500_100000cm-1.tran"
    #fname1 = "PCb_Venus_10bar_TRAN.tran"
    #fname2 = "PCb_Venus_90bar_TRAN.tran"
    fname1 = "PCb_Venus_10bar.trnst"
    fname2 = "PCb_Venus_90bar.trnst"
    title  = "Venus-like (Clouds)"
    title1 = "10 bar"
    title2 = "90 bar"
    plot_kwargs = [
        {"color" : "purple", "label" : title1, "alpha" : 0.7},
        {"color" : "orange", "label" : title2, "alpha" : 0.7}
    ]
    lammin = 0.2
    lammax = 20.
    ylim = [0.0, 100.0]
    legloc = 2
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]
    m = {
        "CO$_2$" : [(2.01, 67.0), (2.7, 77.0), (4.2, 82.0), (5.0, 67.0), (15.0, 82.0)],
        "H$_2$SO$_4$" : [(2.8, 27.), (6.1, 59.), (9.0, 59.)]
    }

    #"H$_2$SO$_4$" : [(6.1, 59.)]

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
                 plot_kwargs=plot_kwargs, legend=True, moleloc=[m, None], legloc=legloc)

    return
#########################################

if __name__ == "__main__":

    make_fig()
