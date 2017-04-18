"""
Figure 28:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig28.5"
    fname1 = "HAZE_1.00e-02ch4_clear_new.trnst"
    fname2 = "HAZE_1.00e-02ch4_cirrus_new.trnst"
    fname3 = "HAZE_1.00e-02ch4_stcum_new.trnst"
    fname4 = "HAZE_1.50e-02ch4_clear_new.trnst"
    fname5 = "HAZE_1.50e-02ch4_cirrus_new.trnst"
    fname6 = "HAZE_1.50e-02ch4_stcum_new.trnst"
    fnames = [fname1, fname2, fname3, fname4, fname5, fname6]
    title  = "Archean Earth-like"
    title1 = "No Haze, Clear"
    title2 = "No Haze, Cirrus Cloud"
    title3 = "No Haze, Stratocumulus Cloud"
    title4 = "Haze, Clear"
    title5 = "Haze, Cirrus Cloud"
    title6 = "Haze, Stratocumulus Cloud"
    plot_kwargs = [
        {"color" : "C1", "label" : title1, "alpha" : 0.7},
        {"color" : "C2", "label" : title2, "alpha" : 0.7},
        {"color" : "C3", "label" : title3, "alpha" : 0.7},
        {"color" : "C4", "label" : title4, "alpha" : 0.7},
        {"color" : "C5", "label" : title5, "alpha" : 0.7},
        {"color" : "C5", "label" : title6, "alpha" : 0.7}
    ]
    lammin = 0.21
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
    fname = [os.path.join(os.path.dirname(__file__),"model_outputs/all/", f) for f in fnames]
    #title = [title1, title2]
    title = [title, "", "", "", "", ""]
    labels = [labels1, labels2]

    # Make reflectance plot
    spc.plot_trans(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, forced_single=True,
                 plot_kwargs=plot_kwargs, legend=True)

    return
#########################################

if __name__ == "__main__":

    make_fig()
