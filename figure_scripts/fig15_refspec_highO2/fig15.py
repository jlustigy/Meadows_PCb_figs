"""
Figure 15:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig15"
    #fname1 = "profile_o2lb_10bar_h2o.pt_filtered_hitran2012_50_100000cm_toa.rad"
    #fname2 = "profile_o2lb_10bar_dry.pt_filtered_hitran2012_50_100000cm_toa.rad"
    fname1 = "10bar_O2_wet.pt_filtered_hitran2012_50_100000cm_toa.rad"
    fname2 = "10bar_O2_dry.pt_filtered_hitran2012_50_100000cm_toa.rad"
    title1 = "10 bar O$_2$ (Ocean)"
    title2 = "10 bar O$_2$ (Desiccated)"
    lammin = 0.15
    lammax = 2.5
    ylim = [-0.02, 0.44]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]

    m1 = {"O$_2$" : [(0.2, 0.1), (0.71, 0.14), (0.77, 0.115)],
          "CO" : [(0.2, 0.05)],
          "O$_3$" : [(0.22, 0.15)],
          "O$_4$" : [(0.39, 0.40),(0.46, 0.34), (0.49, 0.31), (0.54, 0.27), (0.59, 0.22), (0.64, 0.18), (1.06, 0.025), (1.27, 0.025)],
          "H$_2$O" : [(0.835, .09), (0.955, 0.07), (1.4, 0.025), (1.85, 0.025), (2.4, 0.025)],
          "CO$_2$" : [(1.6, 0.06), (2.01, 0.025), (2.2, 0.06)]
          }
    m2 = {"O$_2$" : [(0.2, 0.1), (0.69, 0.23), (0.76, 0.23)],
          "CO" : [(0.2, 0.05), (2.35, 0.15)],
          "O$_3$" : [(0.22, 0.15)],
          "O$_4$" : [(0.39, 0.395), (0.46, 0.34), (0.49, 0.31), (0.54, 0.27), (0.57, 0.15), (0.63, 0.15), (1.06, 0.025), (1.27, 0.025)],
          "CO$_2$" : [(1.6, 0.25), (2.01, 0.1), (2.2, 0.18)]
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
    spc.plot_rad(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, moleloc=moleloc)

    return
#########################################

if __name__ == "__main__":

    make_fig()
