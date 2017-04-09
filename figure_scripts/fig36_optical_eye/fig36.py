"""
Figure 36

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig36_new"
    fname = "fig35_modern_earth_standard_clearsky_2000_100000cm_60sza_sur.rad"
    title  = ""

    modir = "model_outputs/"

    files = [
        {
            "name" : "Proxima Centauri",
            "file" : modir+"proxima_cen_sed.txt",
            "header" : 25,
            "iwlflx" : (0,1)
        },
        {
            "name" : "Earth-like",
            "file" : modir+"profile_Earth_proxb_clear_toa.rad",
            "header" : 0,
            "iwlflx" : (0,3)
        },
        {
            "name" : "Hazy Archean Earth-like",
            "file" : modir+"HAZE_1.50e-02ch4_clear_new_toa.rad",
            "header" : 0,
            "iwlflx" : (0,3)
        },
        {
            "name" : "O$_2$-dominated (Ocean)",
            "file" : modir+"profile_o2lb_10bar_h2o.pt_filtered_hitran2012_50_100000cm_toa.rad",
            "header" : 0,
            "iwlflx" : (0,3)
        },
        {
            "name" : "O$_2$-dominated (Dessicated)",
            "file" : modir+"profile_o2lb_10bar_dry.pt_filtered_hitran2012_50_100000cm_toa.rad",
            "header" : 0,
            "iwlflx" : (0,3)
        },
        {
            "name" : "Venus-like",
            "file" : modir+"PCb_Venus_90bar_toa.rad",
            "header" : 0,
            "iwlflx" : (0,3)
        }
    ]

    lammin = 0.0
    lammax = 10.0
    ylim = [0.0, 1.0]

    # import basic dependencies
    import numpy as np
    import os
    import sys

    # import local package
    sys.path.append("../../")
    from utils import fig_params
    import utils.spectra as spc
    import utils.eyecolors as eye

    # More general params
    #fname = os.path.join(os.path.dirname(__file__),"model_outputs/", fname)

    # Plot and save
    eye.plot_multipanel_ispec(files, savetag=savetag)

    return
#########################################

if __name__ == "__main__":

    make_fig()
