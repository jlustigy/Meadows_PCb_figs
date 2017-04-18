"""
Figure 21:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # import basic dependencies
    import numpy as np
    import os
    import sys

    # import local package
    sys.path.append("../../")
    from utils import fig_params
    import utils.spectra as spc

    # Params specific to this plot
    savetag = "fig21"
    title1 = "Archean Earth-like ($50\%$ clouds)"
    plot_kwargs = [
        {"color" : "purple", "label" : "No Haze, $5\%$ CO$_2$, $1\%$ CH$_4$", "alpha" : 0.7},
        {"color" : "orange", "label" : "Haze, $5\%$ CO$_2$, $1.5\%$ CH$_4$", "alpha" : 0.7}
    ]
    #title2 = "10 bar O$_2$ (Dessicated)"
    lammin = 0.25
    lammax = 2.5
    ylim = [-0.02, 0.3]
    fname_new1 = "model_outputs/HAZE_1.00e-02ch4_combined_new_toa.rad"
    fname_new2 = "model_outputs/HAZE_1.50e-02ch4_combined_new_toa.rad"
    fname_new = [fname_new1, fname_new2]

    m = {
        "Haze" : [(0.35, 0.125)],
        "H$_2$O" : [(0.72, 0.05), (0.815, 0.07), (0.94, 0.19), (1.163, 0.16), (1.47, 0.12), (1.87, 0.06)],
        "CH$_4$" : [(0.62, 0.13), (0.79, 0.2), (0.89, 0.175), (1.0, 0.175), (1.22, 0.19), (1.4, 0.15), (1.7, 0.06), (2.37,0.04)],
        "CO$_2$" : [(1.6, 0.135), (2.05,0.06)]
    }

    fname1 = "HAZE_1.00e-02ch4_clear_new_toa.rad"
    fname2 = "HAZE_1.00e-02ch4_cirrus_new_toa.rad"
    fname3 = "HAZE_1.00e-02ch4_stcum_new_toa.rad"
    fname4 = "HAZE_1.50e-02ch4_clear_new_toa.rad"
    fname5 = "HAZE_1.50e-02ch4_cirrus_new_toa.rad"
    fname6 = "HAZE_1.50e-02ch4_stcum_new_toa.rad"
    fname = [
        os.path.join(os.path.dirname(__file__),"model_outputs/", f)
        for f in [fname1, fname2, fname3, fname4, fname5, fname6]
    ]

    # Weight each flux
    dic1 = {
        fname[0] : 0.5,
        fname[1] : 0.25,
        fname[2] : 0.25
    }
    dic2 = {
        fname[3] : 0.5,
        fname[4] : 0.25,
        fname[5] : 0.25
    }
    # Weight spectra and create output file
    spc.weight_spectra(dic1, output=fname_new1, ret=False)
    spc.weight_spectra(dic2, output=fname_new2, ret=False)

    # More general params
    title = [title1, ""]

    # Make reflectance plot
    spc.plot_rad(fname_new, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, forced_single=True,
                 plot_kwargs=plot_kwargs, legend=True, moleloc=[m, None])

    return
#########################################

if __name__ == "__main__":

    make_fig()
