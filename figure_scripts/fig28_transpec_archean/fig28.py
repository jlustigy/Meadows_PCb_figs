"""
Figure 28:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig28"
    fname1 = "fig28_HAZE_msun21_0.0Ga_1.50e-02ch4_rmix_5.0E-2TRAN.tran"
    fname2 = "fig28_HAZE_msun21_0.0Ga_1.00e-02ch4_rmix_5.0E-2__30.66fscaleTRAN.tran"
    title  = "Archean Earth-like"
    title1 = "$5\%$ CO$_2$, $1.5\%$ CH$_4$ (Haze)"
    title2 = "$5\%$ CO$_2$, $1\%$ CH$_4$ (No Haze)"
    plot_kwargs = [
        {"color" : "orange", "label" : title1, "alpha" : 0.7},
        {"color" : "purple", "label" : title2, "alpha" : 0.7}
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
    fname = [os.path.join(os.path.dirname(__file__),"model_outputs/", f) for f in [fname1, fname2]]
    #title = [title1, title2]
    title = [title, ""]
    labels = [labels1, labels2]

    # Make reflectance plot
    spc.plot_trans(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, forced_single=True,
                 plot_kwargs=plot_kwargs, legend=True)

    return
#########################################

if __name__ == "__main__":

    make_fig()
