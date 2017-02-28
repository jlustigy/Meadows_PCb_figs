"""
Figure 17:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig17"
    fname1 = "fig17_smart_spectra_pandora10bar_cloudy_500_100000cm-1_toa.rad"
    fname2 = "fig17_smart_spectra_pandora90bar_clouds_500_100000cm-1_toa.rad"
    title1 = "10 bar Venus-like (Clouds)"
    title2 = "90 bar Venus-like (Clouds)"
    lammin = 0.3
    lammax = 2.5
    ylim = [-0.02, 0.61]
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
    title = [title1, title2]
    labels = [labels1, labels2]

    # Make reflectance plot
    spc.plot_rad(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None)

    return
#########################################

if __name__ == "__main__":

    make_fig()
