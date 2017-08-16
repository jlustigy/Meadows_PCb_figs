"""
Figure 35

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig35"
    fname = "fig35_modern_earth_standard_clearsky_2000_100000cm_60sza_sur.rad"
    title  = ""

    modir = "model_outputs/"

    files = [
        {
            "name" : "Proxima Centauri (star)",
            "file" : modir+"proxima_cen_sed.txt",
            "header" : 25,
            "iwlflx" : (0,1)
        },
        {
            "name" : "Earth-like (clearsky)",
            "file" : modir+"profile_Earth_proxb_.pt_hitran2012_o4_noh2co_187Kstrat_toa.rad",
            "header" : 0,
            "iwlflx" : (0,3)
        },
        {
            "name" : "Hazy Archean Earth-like (clearsky)",
            "file" : modir+"HAZE_1.50e-02ch4_clear_new_toa.rad",
            "header" : 0,
            "iwlflx" : (0,3)
        },
        {
            "name" : "O$_2$-dominated (ocean)",
            "file" : modir+"10bar_O2_wet.pt_filtered_hitran2012_50_100000cm_toa.rad",
            "header" : 0,
            "iwlflx" : (0,3)
        },
        {
            "name" : "O$_2$-dominated (desiccated)",
            "file" : modir+"10bar_O2_dry.pt_filtered_hitran2012_50_100000cm_toa.rad",
            "header" : 0,
            "iwlflx" : (0,3)
        },
        {
            "name" : "Venus-like (cloudy)",
            "file" : modir+"PCb_Venus_90bar_toa.rad",
            "header" : 0,
            "iwlflx" : (0,3)
        }
    ]

    """
    {
        "name" : "O$_2$-dominated (Ocean)",
        "file" : modir+"10bar_O2_wet.pt_filtered_hitran2012_50_100000cm_toa.rad",
        "header" : 0,
        "iwlflx" : (0,3)
    },
    {
        "name" : "O$_2$-dominated (Dessicated)",
        "file" : modir+"10bar_O2_dry.pt_filtered_hitran2012_50_100000cm_toa.rad",
        "header" : 0,
        "iwlflx" : (0,3)
    },
    {
        "name" : "1 bar CO$_2$/CO/O$_2$ (Dessicated)",
        "file" : modir+"smart_gao_1bar_update_xsec_toa.rad",
        "header" : 0,
        "iwlflx" : (0,3)
    },
    {
        "name" : "10 bar O$_2$-CO$_2$",
        "file" : modir+"10bar_O2_CO2_final.pt_filtered_hitran2012_50_100000cm_toa.rad",
        "header" : 0,
        "iwlflx" : (0,3)
    },
    """

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
